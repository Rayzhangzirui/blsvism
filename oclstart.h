#ifndef OCLSTART_H
#define OCLSTART_H

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include <string>
const cl_uint WORKING_PLATFORM_INDEX = 0;

using namespace std;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void errchk(string thefunction, cl_int errcode_ret);

template <class theType>
string tostring(theType tobeconverted);

template <class theType>
string outputdatatype(theType *root, int numentries);

string outputdatatype(char *root, int numentries);

template <class datatype>
string getdeviceinfo(cl_device_id device, cl_device_info param_name);

string getplatforminfo(cl_platform_id platform, cl_platform_info param_name);

string getallplatforminfo(cl_platform_id platform);

string getsomedeviceinfo(cl_device_id device);

template <class theType>
string getkernerlwginfo(cl_kernel kernel, cl_device_id device, 
                        cl_kernel_work_group_info param_name);

string getallkernelwginfo(cl_kernel kernel, cl_device_id device);

cl_platform_id getPlatform();

cl_device_id *getallDevices(cl_platform_id platform, cl_uint &num_devices);

cl_device_id getDevice(cl_platform_id platform);

cl_device_id getDevice(cl_platform_id platform, cl_device_type deviceType);

cl_device_id getthisDevice(cl_platform_id platform, cl_uint deviceType);


cl_context getContext(cl_device_id device);

cl_command_queue getCommandQueue(cl_context thecontext, cl_device_id thedevice);

cl_program getProgram(const char *kernelFile,
                      cl_context thecontext,
                      cl_device_id thedevice,
                      const char *options);

cl_kernel getKernel(cl_program theprogram, const char *programname);

void getProgramBinaries(cl_program theprogram);


#endif//header guard