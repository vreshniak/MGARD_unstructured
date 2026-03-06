/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * CompressMGARDMeshSerializeOperator.cpp :
 *
 *  Created on: Dec 1, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "CompressMGARDMeshSerializeOperator.h"

#include "adios2.h"
#include "adios2/helper/adiosFunctions.h"
#include <mgard/MGARDConfig.hpp>
#include <mgard/compress_x.hpp>

#include <cstring>
#include <mutex>


template <typename T>
auto reorder(const std::vector<T> &vec, const std::vector<int64_t> &order)
-> std::vector<T>
{
    std::vector<T> new_vec(vec.size());

    for (size_t i=0; i<vec.size(); i++)
        new_vec[i] = vec[order[i]];

    return new_vec;
}


template <typename T>
auto inverse_reorder(const std::vector<T> &vec, const std::vector<int64_t> &order)
-> std::vector<T>
{
    std::vector<T> new_vec(vec.size());

    for (size_t i=0; i<vec.size(); i++)
        new_vec[order[i]] = vec[i];

    return new_vec;
}



std::vector<int64_t> node_order;
size_t curr_bId = -1;

bool verbose = false;
bool debugging = false;


namespace adios2
{
namespace plugin
{

/* "STATIC" PART of operator to have read mesh only once */
std::string meshFileName;
std::mutex readMeshMutex;
bool meshReadSuccessfully = false;
// define mesh mapping variable files


void ReadMesh(std::string meshfile, std::string node_order_name, size_t bId) //, std::string connectivity_name, std::vector<std::string> coo_names, size_t bId, int n_nodes_in_element) //, std::string Blc)
{
    meshReadSuccessfully = true;
    meshFileName = meshfile;

    if (debugging) std::cout << "Debugging: Reading Mesh File " << meshfile << std::endl; // <<  "starting from block " << bId << "\n";

    adios2::ADIOS adios;
    adios2::IO reader_io = adios.DeclareIO("BPReader");
    reader_io.SetEngine("BPFile");
    adios2::Engine bpReader = reader_io.Open(meshFileName, adios2::Mode::Read);

    // size_t offsetGrid = 0, offsetNode = 0;
    while (true) {
        adios2::StepStatus read_status = bpReader.BeginStep(adios2::StepMode::Read, -1.0f);
        if (read_status == adios2::StepStatus::NotReady) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            continue;
        }
        else if (read_status != adios2::StepStatus::OK) {
            break;
        }

        // read node ordering from meshfile
        // std::vector<int64_t> node_order;
        node_order.clear();
        adios2::Variable<int64_t> adios_node_order = reader_io.InquireVariable<int64_t>(node_order_name);
        adios_node_order.SetBlockSelection(bId);
        bpReader.Get(adios_node_order, node_order, adios2::Mode::Sync);

        // std::vector<int64_t> node_order;
        // adios2::Variable<int64_t> adios_node_order = reader_io.InquireVariable<int64_t>("/hpMusic_base/hpMusic_Zone/Elem/ElementConnectivity");
        // adios_node_order.SetBlockSelection(curr_bId);
        // bpReader.Get(adios_node_order, node_order, adios2::Mode::Sync);
        // break;


        // // allocate memory for reading variables
        // std::vector<int64_t> ElementConnectivity;
        // std::vector<std::vector<double>> coos(coo_names.size());

        // adios2::Variable<int64_t> adios_connectivity = reader_io.InquireVariable<int64_t>(connectivity_name);
        // std::vector<adios2::Variable<double>> adios_coos(coo_names.size());
        // for (int i=0; i<coo_names.size(); i++)
        //     adios_coos[i] = reader_io.InquireVariable<double>(coo_names[i]);

        // // set block/subregion for the variables
        // adios_connectivity.SetBlockSelection(bId);
        // for (auto &coo : adios_coos) coo.SetBlockSelection(bId);

        // // read variables
        // bpReader.Get(adios_connectivity, ElementConnectivity, adios2::Mode::Sync);
        // for (int i=0; i<coo_names.size(); i++) bpReader.Get(adios_coos[i], coos[i], adios2::Mode::Sync);

        // bpReader.PerformGets();

        // //////////////////////////////////////////////////////////////////////////////////
        // // reordering

        // size_t num_nodes = coos[0].size();

        // // compute connectivity of the node graph
        // auto GlobalConnectivity = ComputeGlobalConnectivity(ElementConnectivity, num_nodes, n_nodes_in_element);

        // node_order = find_node_order(GlobalConnectivity, coos);

        // std::vector<size_t> reorder_map(node_order.size());
        // for (size_t i=0; i<reorder_map.size(); i++)
        //     reorder_map[node_order[i]] = i;

        // for (size_t i=0; i<ElementConnectivity.size(); i++)
        //     ElementConnectivity[i] = reorder_map[ElementConnectivity[i]];

        bpReader.EndStep();
    }
    bpReader.Close();
    return;
}

/* END OF "STATIC" PART */


CompressMGARDMeshSerializeOperator::CompressMGARDMeshSerializeOperator(const Params &parameters) : PluginOperatorInterface(parameters)
{
    if (verbose) std::cout << "=== CompressMGARDMeshSerializeOperator constructor ===" << std::endl;
    auto itMeshFileName  = m_Parameters.find("meshfile");
    auto itNodeOrderName = m_Parameters.find("node_order");
    if (itMeshFileName == m_Parameters.end())
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressMGARDMeshSerializeOperator", "constructor", "This operator needs an unstructured mesh input file with parameter name 'meshfile'");
    }
    std::string &meshfile = itMeshFileName->second;

    auto itBlockId = m_Parameters.find("blockid");
    if (itBlockId == m_Parameters.end())
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressMGARDMeshSerializeOperator", "constructor", "This operator needs a blockId input 'blockid'");
    }
    size_t bId = (size_t)std::stoi(itBlockId->second);
    if (!meshReadSuccessfully || curr_bId!=bId)
    {
        std::lock_guard<std::mutex> lck(readMeshMutex);
        if (!meshFileName.empty() && meshfile != meshFileName)
        {
            helper::Throw<std::invalid_argument>("Operator", "CompressMGARDMeshSerializeOperator", "constructor", "Cannot process more than one mesh files. Already read " + meshFileName);
        }
        curr_bId = bId;
        ReadMesh(itMeshFileName->second, itNodeOrderName->second, bId); //, "connectivity", {"x","y","z"}, bId, 8);
    }
}


size_t CompressMGARDMeshSerializeOperator::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount, const DataType type, char *bufferOut)
{
    if (verbose) std::cout << "=== CompressMGARDMeshSerializeOperator::Operate() ===" << std::endl;

    Dims convertedDims = ConvertDims(blockCount, type, 3);
    const size_t ndims = convertedDims.size();
    if (ndims > 5) helper::Throw<std::invalid_argument>("Operator", "CompressMGARDMeshSerializeOperator", "Operate", "MGARD does not support data in " + std::to_string(ndims) + " dimensions");

    ////////////////////////////////////////////////////////////////////////////////
    // mgard V1 metadata

    const uint8_t bufferVersion = 1;
    size_t bufferOutOffset = 0;

    MakeCommonHeader(bufferOut, bufferOutOffset, bufferVersion);

    PutParameter(bufferOut, bufferOutOffset, ndims);
    for (const auto &d : convertedDims)
        PutParameter(bufferOut, bufferOutOffset, d);
    PutParameter(bufferOut, bufferOutOffset, type);
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_MAJOR));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_MINOR));
    PutParameter(bufferOut, bufferOutOffset, static_cast<uint8_t>(MGARD_VERSION_PATCH));

    // mgard V1 metadata end
    ////////////////////////////////////////////////////////////////////////////////


    // set MGARD type
    mgard_x::data_type mgardType;
    if (type == helper::GetDataType<float>())
        mgardType = mgard_x::data_type::Float;
    else if (type == helper::GetDataType<double>())
        mgardType = mgard_x::data_type::Double;
    else if (type == helper::GetDataType<std::complex<float>>())
        mgardType = mgard_x::data_type::Float;
    else if (type == helper::GetDataType<std::complex<double>>())
        mgardType = mgard_x::data_type::Double;
    else
        helper::Throw<std::invalid_argument>("Operator", "CompressMGARDMeshSerializeOperator", "Operate", "MGARD only supports float and double types");
    // set MGARD type end


    // set MGARD style dim info
    mgard_x::DIM mgardDim = ndims;
    std::vector<mgard_x::SIZE> mgardCount;
    for (const auto &c : convertedDims)
        mgardCount.push_back(c);
    // set MGARD style dim info end

    // Parameters
    bool hasTolerance = false;
    double tolerance = 0.0;
    double s = 0.0;
    auto errorBoundType = mgard_x::error_bound_type::REL;

    // input size under this bound will not compress
    size_t thresholdSize = 100000;

    auto itThreshold = m_Parameters.find("threshold");
    if (itThreshold != m_Parameters.end())
    {
        thresholdSize = std::stod(itThreshold->second);
    }
    auto itAccuracy = m_Parameters.find("accuracy");
    if (itAccuracy != m_Parameters.end())
    {
        tolerance = std::stod(itAccuracy->second);
        hasTolerance = true;
    }
    auto itTolerance = m_Parameters.find("tolerance");
    if (itTolerance != m_Parameters.end())
    {
        tolerance = std::stod(itTolerance->second);
        hasTolerance = true;
    }
    if (!hasTolerance)
    {
        helper::Throw<std::invalid_argument>("Operator", "CompressMGARDMeshSerializeOperator", "Operate", "missing mandatory parameter tolerance / accuracy");
    }
    auto itSParameter = m_Parameters.find("s");
    if (itSParameter != m_Parameters.end())
    {
        s = std::stod(itSParameter->second);
    }
    auto itMode = m_Parameters.find("mode");
    if (itMode != m_Parameters.end())
    {
        if (itMode->second == "ABS")
        {
            errorBoundType = mgard_x::error_bound_type::ABS;
        }
        else if (itMode->second == "REL")
        {
            errorBoundType = mgard_x::error_bound_type::REL;
        }
    }

    // let mgard know the output buffer size
    size_t sizeOut = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type));

    if (sizeOut < thresholdSize)
    {
        /* disable compression and add marker in the header*/
        PutParameter(bufferOut, bufferOutOffset, false);
        headerSize = bufferOutOffset;
        return 0;
    }
    else
    {
        /* enable compression and add marker in the header*/
        PutParameter(bufferOut, bufferOutOffset, true);
    }


    ///////////////////////////////////////////////////////
    // compress

    mgard_x::Config config;
    config.lossless = mgard_x::lossless_type::Huffman_Zstd;

    void *compressedData = bufferOut + bufferOutOffset;
    // mgard_x::compress(mgardDim, mgardType, mgardCount, tolerance, s, errorBoundType, dataIn, compressedData, sizeOut, config, true);
    std::vector<double> datvec((double*)dataIn, (double*)dataIn+mgardCount[0]);
    mgard_x::compress(mgardDim, mgardType, mgardCount, tolerance, s, errorBoundType, reorder(datvec,node_order).data(), compressedData, sizeOut, config, true);
    bufferOutOffset += sizeOut;

    // end compress
    ///////////////////////////////////////////////////////

    return bufferOutOffset;
}

size_t CompressMGARDMeshSerializeOperator::GetHeaderSize() const { return headerSize; }


size_t CompressMGARDMeshSerializeOperator::DecompressV1(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    // Do NOT remove even if the buffer version is updated. Data might be still
    // in lagacy formats. This function must be kept for backward compatibility.
    // If a newer buffer format is implemented, create another function, e.g.
    // DecompressV2 and keep this function for decompressing lagacy data.

    size_t bufferInOffset = 0;

    ////////////////////////////////////////////////////////////////////////////////
    // read header
    const size_t ndims = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    Dims blockCount(ndims);
    for (size_t i = 0; i < ndims; ++i)
        blockCount[i] = GetParameter<size_t, size_t>(bufferIn, bufferInOffset);
    const DataType type = GetParameter<DataType>(bufferIn, bufferInOffset);
    m_VersionInfo = " Data is compressed using MGARD Version " + std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                    "." + std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) + "." +
                    std::to_string(GetParameter<uint8_t>(bufferIn, bufferInOffset)) +
                    ". Please make sure a compatible version is used for decompression.";

    const bool isCompressed = GetParameter<bool>(bufferIn, bufferInOffset);
    // read header end
    ////////////////////////////////////////////////////////////////////////////////


    // decompressed size
    size_t sizeOut = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type));
    if (type == DataType::FloatComplex || type == DataType::DoubleComplex)
        sizeOut /= 2;

    if (isCompressed)
    {
        try
        {
            void *dataOutVoid = dataOut;
            mgard_x::decompress(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOutVoid, true);

            if (type == helper::GetDataType<float>()){
                std::vector<float> datvec((float*)dataOutVoid, (float*)dataOutVoid+sizeOut/sizeof(float));
                dataOutVoid = (void*)inverse_reorder(datvec,node_order).data();
            }
            else if (type == helper::GetDataType<double>()){
                std::vector<double> datvec((double*)dataOutVoid, (double*)dataOutVoid+sizeOut/sizeof(double));
                memcpy(dataOutVoid, (void*)inverse_reorder(datvec,node_order).data(), sizeOut);
                // dataOutVoid = (void*)inverse_reorder(datvec,node_order).data();
            }
            // std::vector<double> datvec((double*)dataIn, (double*)dataIn+sizeIn-bufferInOffset);
        }
        catch (...)
        {
            helper::Throw<std::runtime_error>("Operator", "CompressMGARDMeshSerializeOperator", "DecompressV1", m_VersionInfo);
        }
        return sizeOut;
    }

    headerSize += bufferInOffset;
    return 0;
}

size_t CompressMGARDMeshSerializeOperator::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    size_t bufferInOffset = 1; // skip operator type
    const uint8_t bufferVersion = GetParameter<uint8_t>(bufferIn, bufferInOffset);
    bufferInOffset += 2; // skip two reserved bytes
    headerSize = bufferInOffset;

    if (bufferVersion == 1)
    {
        return DecompressV1(bufferIn + bufferInOffset, sizeIn - bufferInOffset, dataOut);
    }
    else if (bufferVersion == 2)
    {
        // TODO: if a Version 2 mgard buffer is being implemented, put it here
        // and keep the DecompressV1 routine for backward compatibility
    }
    else
    {
        helper::Throw<std::runtime_error>("Operator", "CompressMGARDMeshSerializeOperator", "InverseOperate", "invalid mgard buffer version");
    }

    return 0;
}

bool CompressMGARDMeshSerializeOperator::IsDataTypeValid(const DataType type) const
{
    if (type == DataType::Double || type == DataType::Float || type == DataType::DoubleComplex || type == DataType::FloatComplex)
    {
        return true;
    }
    return false;
}

} // end namespace plugin
} // end namespace adios2

extern "C" {

adios2::plugin::CompressMGARDMeshSerializeOperator *OperatorCreate(const adios2::Params &parameters)
{
    return new adios2::plugin::CompressMGARDMeshSerializeOperator(parameters);
}

void OperatorDestroy(adios2::plugin::CompressMGARDMeshSerializeOperator *obj) { delete obj; }
}
