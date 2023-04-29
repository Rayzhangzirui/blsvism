#define CL_SILENCE_DEPRECATION //opencl deprecated in macOS 10.14 Mojave
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>

#include "oclstart.h"

////////////////////////////////////////////////////////////////////////////

using namespace std;

//error codes
const char *err_code (cl_int err_in)
{
    switch (err_in) {
        case CL_SUCCESS:
            return (char*)"CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND:
            return (char*)"CL_DEVICE_NOT_FOUND";
        case CL_DEVICE_NOT_AVAILABLE:
            return (char*)"CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE:
            return (char*)"CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            return (char*)"CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES:
            return (char*)"CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY:
            return (char*)"CL_OUT_OF_HOST_MEMORY";
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            return (char*)"CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP:
            return (char*)"CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH:
            return (char*)"CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            return (char*)"CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE:
            return (char*)"CL_BUILD_PROGRAM_FAILURE";
        case CL_MAP_FAILURE:
            return (char*)"CL_MAP_FAILURE";
        case CL_MISALIGNED_SUB_BUFFER_OFFSET:
            return (char*)"CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
            return (char*)"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        case CL_INVALID_VALUE:
            return (char*)"CL_INVALID_VALUE";
        case CL_INVALID_DEVICE_TYPE:
            return (char*)"CL_INVALID_DEVICE_TYPE";
        case CL_INVALID_PLATFORM:
            return (char*)"CL_INVALID_PLATFORM";
        case CL_INVALID_DEVICE:
            return (char*)"CL_INVALID_DEVICE";
        case CL_INVALID_CONTEXT:
            return (char*)"CL_INVALID_CONTEXT";
        case CL_INVALID_QUEUE_PROPERTIES:
            return (char*)"CL_INVALID_QUEUE_PROPERTIES";
        case CL_INVALID_COMMAND_QUEUE:
            return (char*)"CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_HOST_PTR:
            return (char*)"CL_INVALID_HOST_PTR";
        case CL_INVALID_MEM_OBJECT:
            return (char*)"CL_INVALID_MEM_OBJECT";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            return (char*)"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case CL_INVALID_IMAGE_SIZE:
            return (char*)"CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_SAMPLER:
            return (char*)"CL_INVALID_SAMPLER";
        case CL_INVALID_BINARY:
            return (char*)"CL_INVALID_BINARY";
        case CL_INVALID_BUILD_OPTIONS:
            return (char*)"CL_INVALID_BUILD_OPTIONS";
        case CL_INVALID_PROGRAM:
            return (char*)"CL_INVALID_PROGRAM";
        case CL_INVALID_PROGRAM_EXECUTABLE:
            return (char*)"CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_KERNEL_NAME:
            return (char*)"CL_INVALID_KERNEL_NAME";
        case CL_INVALID_KERNEL_DEFINITION:
            return (char*)"CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL:
            return (char*)"CL_INVALID_KERNEL";
        case CL_INVALID_ARG_INDEX:
            return (char*)"CL_INVALID_ARG_INDEX";
        case CL_INVALID_ARG_VALUE:
            return (char*)"CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE:
            return (char*)"CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS:
            return (char*)"CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_WORK_DIMENSION:
            return (char*)"CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE:
            return (char*)"CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE:
            return (char*)"CL_INVALID_WORK_ITEM_SIZE";
        case CL_INVALID_GLOBAL_OFFSET:
            return (char*)"CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_EVENT_WAIT_LIST:
            return (char*)"CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_EVENT:
            return (char*)"CL_INVALID_EVENT";
        case CL_INVALID_OPERATION:
            return (char*)"CL_INVALID_OPERATION";
        case CL_INVALID_GL_OBJECT:
            return (char*)"CL_INVALID_GL_OBJECT";
        case CL_INVALID_BUFFER_SIZE:
            return (char*)"CL_INVALID_BUFFER_SIZE";
        case CL_INVALID_MIP_LEVEL:
            return (char*)"CL_INVALID_MIP_LEVEL";
        case CL_INVALID_GLOBAL_WORK_SIZE:
            return (char*)"CL_INVALID_GLOBAL_WORK_SIZE";
        case CL_INVALID_PROPERTY:
            return (char*)"CL_INVALID_PROPERTY";

        default:
            return (char*)"UNKNOWN ERROR";
    }
}
////////////////////////////////////////////////////////////////////////////

void errchk(string thefunction, cl_int errcode_ret)
{
   if (errcode_ret != CL_SUCCESS)
   {
      cout << thefunction << " returned error: " << err_code(errcode_ret)<< "\n";
      exit(1);
   }
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class theType>
string tostring(theType tobeconverted)
{
   stringstream ss;
   ss << tobeconverted;
   return ss.str();
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class theType>
string outputdatatype(theType *root, int numentries)
{
   string ans;

   ans = tostring<theType>(root[0]);

   for (int i=1; i<numentries; ++i)
   {
      ans += " ";
      ans += tostring<theType>(root[i]);
   }

   return ans; 

}

string outputdatatype(char *root, int numentries)
{
   string ans;

   for (int i=0; i<numentries; ++i)
      ans += tostring(root[i]);

   return ans;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class datatype>
string getdeviceinfo(cl_device_id device, cl_device_info param_name)
{
   size_t param_value_size_ret;
   clGetDeviceInfo(device, param_name, 0, NULL, &param_value_size_ret);

   size_t numentries = param_value_size_ret/sizeof(datatype);

   datatype *param_value = new datatype[numentries];

   clGetDeviceInfo(device, param_name, param_value_size_ret, param_value, NULL);

   string ans;

  ans += outputdatatype(param_value, numentries);

   delete [] param_value;

   return ans;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

string getplatforminfo(cl_platform_id platform, cl_platform_info param_name)
{
   size_t param_value_size_ret = 0;

   // get size of return data:
   clGetPlatformInfo(platform, param_name, 0, NULL, &param_value_size_ret);

   // allocate enough space for return data:
   char *param_value = new char[param_value_size_ret/sizeof(char)];

   // get return data:
   clGetPlatformInfo(platform, param_name, param_value_size_ret,
                     param_value, NULL);

   // convert return data to string:
   string ans = "";
   for (int i=0; i<param_value_size_ret; ++i)
      ans += param_value[i];

   delete [] param_value;

   return ans;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

string getallplatforminfo(cl_platform_id platform)
{
   string ans = "";

   ans += "Platform id: " + tostring<cl_platform_id>(platform);
   ans += "\n";

   ans += "CL_PLATFORM_PROFILE: ";
   ans += getplatforminfo(platform, CL_PLATFORM_PROFILE);
   ans += "\n";

   ans += "CL_PLATFORM_VERSION: ";
   ans += getplatforminfo(platform, CL_PLATFORM_VERSION);
   ans += "\n";

   ans += "CL_PLATFORM_NAME: ";
   ans += getplatforminfo(platform, CL_PLATFORM_NAME);
   ans += "\n";

   ans += "CL_PLATFORM_VENDOR: ";
   ans += getplatforminfo(platform, CL_PLATFORM_VENDOR);
   ans += "\n";

   ans += "CL_PLATFORM_EXTENSIONS: ";
   ans += getplatforminfo(platform, CL_PLATFORM_EXTENSIONS);
   ans += "\n";

   return ans;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

string getsomedeviceinfo(cl_device_id device)
{
   string ans;

   // CL_DEVICE_AVAILABLE:
   ans += "CL_DEVICE_AVAILABLE: ";
   ans += getdeviceinfo<cl_bool>(device, CL_DEVICE_AVAILABLE);
   ans += "\n";

   // CL_DEVICE_COMPILER_AVAILABLE:
   ans += "CL_DEVICE_COMPILER_AVAILABLE: ";
   ans += getdeviceinfo<cl_bool>(device, CL_DEVICE_COMPILER_AVAILABLE);
   ans += "\n";

   // CL_DEVICE_EXTENSIONS:
   ans += "CL_DEVICE_EXTENSIONS: ";
   ans += getdeviceinfo<char>(device, CL_DEVICE_EXTENSIONS);
   ans += "\n";

   // CL_DEVICE_GLOBAL_MEM_CACHE_SIZE:
   ans += "CL_DEVICE_GLOBAL_MEM_CACHE_SIZE: ";
   ans += getdeviceinfo<cl_ulong>(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE);
   ans += "\n";

   // CL_DEVICE_GLOBAL_MEM_SIZE: 
   ans += "CL_DEVICE_GLOBAL_MEM_SIZE: ";
   ans += getdeviceinfo<cl_ulong>(device, CL_DEVICE_GLOBAL_MEM_SIZE);
   ans += "\n";

   // CL_DEVICE_LOCAL_MEM_SIZE: 
   ans += "CL_DEVICE_LOCAL_MEM_SIZE: ";
   ans += getdeviceinfo<cl_ulong>(device, CL_DEVICE_LOCAL_MEM_SIZE);
   ans += "\n";

   // CL_DEVICE_LOCAL_MEM_TYPE: 
   ans += "CL_DEVICE_LOCAL_MEM_TYPE: ";
   ans += getdeviceinfo<cl_device_local_mem_type>(device, 
                                                  CL_DEVICE_LOCAL_MEM_TYPE);
   ans += " (CL_LOCAL: " + tostring<cl_device_local_mem_type>(CL_LOCAL)
        + " CL_GLOBAL: " + tostring<cl_device_local_mem_type>(CL_GLOBAL) + ")";
   ans += "\n";

   // CL_DEVICE_MAX_COMPUTE_UNITS: 
   ans += "CL_DEVICE_MAX_COMPUTE_UNITS: ";
   ans += getdeviceinfo<cl_uint>(device, CL_DEVICE_MAX_COMPUTE_UNITS);
   ans += "\n";

   // CL_DEVICE_MAX_CONSTANT_ARGS:
   ans += "CL_DEVICE_MAX_CONSTANT_ARGS: ";
   ans += getdeviceinfo<cl_uint>(device, CL_DEVICE_MAX_CONSTANT_ARGS);
   ans += "\n";

   // CL_DEVICE_MAX_MEM_ALLOC_SIZE: 
   ans += "CL_DEVICE_MAX_MEM_ALLOC_SIZE: ";
   ans += getdeviceinfo<cl_ulong>(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE);
   ans += "\n";

   // CL_DEVICE_MAX_PARAMETER_SIZE:
   ans += "CL_DEVICE_MAX_PARAMETER_SIZE: ";
   ans += getdeviceinfo<cl_ulong>(device, CL_DEVICE_MAX_PARAMETER_SIZE);
   ans += "\n";

   // CL_DEVICE_MAX_WORK_GROUP_SIZE: 
   ans += "CL_DEVICE_MAX_WORK_GROUP_SIZE: ";
   ans += getdeviceinfo<size_t>(device, CL_DEVICE_MAX_WORK_GROUP_SIZE);
   ans += "\n";

   // CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: 
   ans += "CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS: ";
   ans += getdeviceinfo<cl_uint>(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS);
   ans += "\n";

   // CL_DEVICE_MAX_WORK_ITEM_SIZES: 
   ans += "CL_DEVICE_MAX_WORK_ITEM_SIZES: ";
   ans += getdeviceinfo<size_t>(device, CL_DEVICE_MAX_WORK_ITEM_SIZES);
   ans += "\n";

   // CL_DEVICE_NAME: 
   ans += "CL_DEVICE_NAME: ";
   ans += getdeviceinfo<char>(device, CL_DEVICE_NAME);
   ans += "\n";

   // CL_DEVICE_PLATFORM: 
   ans += "CL_DEVICE_PLATFORM: ";
   ans += getdeviceinfo<cl_platform_id>(device, CL_DEVICE_PLATFORM);
   ans += "\n";

  // CL_DEVICE_PREFFERED_VECTOR_WIDTH_DOUBLE:
   ans += "CL_DEVICE_PREFFERED_VECTOR_WIDTH_DOUBLE: ";
   ans += getdeviceinfo<cl_uint>(device,
                                 CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE);
   ans += "\n";

   // CL_DEVICE_PROFILE:    
   ans += "CL_DEVICE_PROFILE: ";
   ans += getdeviceinfo<char>(device, CL_DEVICE_PROFILE);
   ans += "\n";

   // CL_DEVICE_TYPE: 
   ans += "CL_DEVICE_TYPE: ";
   ans += getdeviceinfo<cl_device_type>(device, CL_DEVICE_TYPE);
   ans += " (device types: CPU: "
          +           tostring<cl_device_type>(CL_DEVICE_TYPE_CPU) + "; "
          + "GPU: " + tostring<cl_device_type>(CL_DEVICE_TYPE_GPU) + "; "
          + "ACCELERATOR: " + tostring<cl_device_type>(CL_DEVICE_TYPE_ACCELERATOR) + "; "
          + "DEFAULT: " + tostring<cl_device_type>(CL_DEVICE_TYPE_DEFAULT) + ")";
   ans += "\n";
                            

   // CL_DEVICE_VENDOR: 
   ans += "CL_DEVICE_VENDOR: ";
   ans += getdeviceinfo<char>(device, CL_DEVICE_VENDOR);
   ans += "\n";

   // CL_DEVICE_VERSION: 
   ans += "CL_DEVICE_VERSION: ";
   ans += getdeviceinfo<char>(device, CL_DEVICE_VERSION);
   ans += "\n";

   // CL_DRIVER_VERSION: 
   ans += "CL_DRIVER_VERSION: ";
   ans += getdeviceinfo<char>(device, CL_DRIVER_VERSION);
   ans += "\n";

   return ans;
  
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class theType>
string getkernelwginfo(cl_kernel kernel, cl_device_id device,
                       cl_kernel_work_group_info param_name)
{
   cl_int err = 0;
   size_t ret_size = 0;

   err = clGetKernelWorkGroupInfo(kernel, device, param_name, 
                                  0, NULL, &ret_size);
   errchk("clGetKernelWorkGroupInfo", err);

   theType *param_value = new theType[ret_size/sizeof(theType)];

   err = clGetKernelWorkGroupInfo(kernel, device, param_name, 
                                  ret_size, param_value, NULL);
   errchk("clGetKernelWorkGroupInfo", err);

   return outputdatatype(param_value, ret_size/sizeof(theType));
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

string getallkernelwginfo(cl_kernel kernel, cl_device_id device)
{
   string ans;

/*
   ans += "CL_KERNEL_GLOBAL_WORK_SIZE: ";
   ans += getkernelwginfo<size_t>(kernel, device, CL_KERNEL_GLOBAL_WORK_SIZE);
   ans += "\n";
*/

   ans += "CL_KERNEL_WORK_GROUP_SIZE: ";
   ans += getkernelwginfo<size_t>(kernel, device, CL_KERNEL_WORK_GROUP_SIZE);
   ans += "\n";

   ans += "CL_KERNEL_COMPILE_WORK_GROUP_SIZE: ";
   ans += getkernelwginfo<size_t>(kernel, device, CL_KERNEL_COMPILE_WORK_GROUP_SIZE);
   ans += "\n";

   ans += "CL_KERNEL_LOCAL_MEM_SIZE: ";
   ans += getkernelwginfo<size_t>(kernel, device, CL_KERNEL_LOCAL_MEM_SIZE);
   ans += "\n";

//   ans += "CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE: ";
//   ans += getkernelwginfo<cl_ulong>(kernel, device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE);
//   ans += "\n";

//   ans += "CL_KERNEL_PRIVATE_MEM_SIZE: ";
//   ans += getkernelwginfo<size_t>(kernel, device, CL_KERNEL_PRIVATE_MEM_SIZE);
//   ans += "\n";

   return ans;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

cl_platform_id getPlatform()
{
   cl_uint num_platforms = 0;

   // get number of platforms available:
   clGetPlatformIDs(0, NULL, &num_platforms);

   //printf("Number of platforms available: %d\n", num_platforms);
   if (num_platforms == 0)
      exit(1);

   // get the platforms:
   cl_platform_id *platforms = new cl_platform_id[num_platforms]; 
   clGetPlatformIDs(num_platforms, platforms, NULL);

   // output platform info:
   // for (cl_uint i=0; i<num_platforms; ++i)
   //    cout << getallplatforminfo(platforms[i]);

   ////////////////////////////////////////////////////////////////////////

   // we work with WORKING_PLATFORM_INDEX:
   cl_platform_id platform;
   if (WORKING_PLATFORM_INDEX < num_platforms)
      platform = platforms[WORKING_PLATFORM_INDEX];
   else
   {
      printf("WORKING_PLATFORM_INDEX = %d is less than num_platforms = %d\n",
              WORKING_PLATFORM_INDEX, num_platforms);
      exit(1);
   }

   ////////////////////////////////////////////////////////////////////////

   delete [] platforms;

   return platform;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

cl_device_id *getallDevices(cl_platform_id platform, cl_uint &num_devices)
{
   // get device id's:

   cl_int err;

   // use CPU only, 
   // use all leads to error <program source>:329:32: error: call to '__fast_relax_log' is ambiguous
   err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices);
   errchk("clGetDeviceIDs", err);

   // printf("\n");
   // printf("Number of devices available: %d\n", num_devices);

   if (num_devices == 0)
      exit(1);

   // get the devices:

   cl_device_id *devices = new cl_device_id[num_devices];
   err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, num_devices, devices, NULL);
   errchk("clGetDeviceIDs", err);
   // output device info
   // for (cl_uint i=0; i<num_devices; ++i)
   // {
   //    printf("\nDevice %d:\n", i);
   //    cout << getsomedeviceinfo(devices[i]);
   // }

   return devices;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/*
cl_device_id getDevice(cl_platform_id platform, cl_uint device_index)
{
   cl_device_id *devices;
   cl_uint num_devices;

   devices = getallDevices(platform, num_devices);

   if (num_devices <= device_index)
   {
      cout << "Requested device index: " << device_index <<
              "larger than number of devices: " << num_devices << ".\n";
      exit(1);
   }

   cout << "\nUsing device number: " << device_index << "\n";

   cl_device_id device;
   device = devices[device_index];
   delete [] devices;

   return device;
}
*/

////////////////////////////////////////////////////////////////////////////
//Not used


// cl_device_id getDevice(cl_platform_id platform)
// {
//    cl_device_id *devices;
//    cl_uint num_devices;

//    devices = getallDevices(platform, num_devices);

//    cout << "\nUsing device number: 0\n";

//    cl_device_id device;
//    device = devices[0];
//    delete [] devices;

//    return device;
// }

////////////////////////////////////////////////////////////////////////////
//Not used

cl_device_id getDevice(cl_platform_id platform, cl_device_type deviceType)
{
   cl_device_id *devices;
   cl_uint num_devices;
   cl_int err;

   devices = getallDevices(platform, num_devices);

   cl_uint i = 0;
   cl_device_type thetype;

   do
   {
      err = clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, sizeof(cl_device_type),
                    &thetype, NULL);
      errchk("clGetDeviceIDs", err);
   }
   while (thetype != deviceType && ++i < num_devices);

   if (i == num_devices)  // no match for desired device type...
   {
      cout << "Error, no device of type " << deviceType << "\n";
      exit(1);
   }

   #ifdef DEBUG
   cout << "\nUsing device number: " << i << "\n";
   #endif

   cl_device_id device;
   device = devices[i];
   delete [] devices;

   return device;
}

//////////////////////////////////////////////////////////////////////

cl_device_id getthisDevice(cl_platform_id platform, cl_uint i)
{
   cl_device_id *devices;
   cl_uint num_devices;

   devices = getallDevices(platform, num_devices);

   #ifdef DEBUG
   cout << "\nUsing device number: " << i << "\n";
   #endif
   cl_device_id device;
   device = devices[i];
   delete [] devices;

   return device;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

cl_context getContext(cl_device_id device)
{
   cl_int errcode_ret;
   cl_context context;

   context = clCreateContext(0, 1, &device, NULL, NULL, &errcode_ret);
   errchk("clCreateContext", errcode_ret);  // check err message

   return context;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

cl_command_queue getCommandQueue(cl_context thecontext, cl_device_id thedevice)
{
   int err;
   cl_command_queue thequeue;

//   thequeue = clCreateCommandQueue(thecontext, thedevice, 0, &err);
   thequeue = clCreateCommandQueue(thecontext, thedevice,
              CL_QUEUE_PROFILING_ENABLE, &err);
   errchk("clCreateCommandQueue", err);

   return thequeue;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

cl_program getProgram(const char *kernelFile,
                      cl_context thecontext,
                      cl_device_id thedevice,
                      const char *options)
{
   int err;
   cl_program theprogram;

   ///////////////////////////////////////////////////////////////////////
   ifstream file(kernelFile);

   string prog(istreambuf_iterator<char>(file), (istreambuf_iterator<char>()));

   file.close();

   ///////////////////////////////////////////////////////////////////////

   char *kersource;

   kersource = new char[prog.size()+1];

   strcpy(kersource, prog.c_str());
   theprogram = clCreateProgramWithSource(thecontext,
                                          1,
                                          (const char**) &kersource,
                                          NULL,
                                          &err);
   errchk("clCreateProgramWithSource", err);

   ///////////////////////////////////////////////////////////////////////

   
   printf("clBuildProgram settings for: ");
   printf("%10.100s: ", kernelFile);
   

   cout << ( (!options) ? "NULL" : options) << endl;
   err = clBuildProgram(theprogram, 1, &thedevice, options, NULL, NULL); 

   if (err != CL_SUCCESS)
   {
      char *buffer;
      size_t len;

      printf("Error (%d) in clBuildProgram\n", err);

      clGetProgramBuildInfo(theprogram, thedevice, CL_PROGRAM_BUILD_LOG,
                            0, buffer, &len);

      buffer = new char [len];

      clGetProgramBuildInfo(theprogram, thedevice, CL_PROGRAM_BUILD_LOG,
                            len, buffer, &len);

      printf("%s\n", buffer);

      delete [] buffer;

      exit(1);
   }

   ///////////////////////////////////////////////////////////////////////

   delete [] kersource;

   return theprogram;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

cl_kernel getKernel(cl_program theprogram, const char *programname)
{
   int err;
   cl_kernel thekernel;

   thekernel = clCreateKernel(theprogram, programname, &err);
   if(err !=CL_SUCCESS){
    printf("failed to build program %10.100s: ", programname); 
   }
   errchk("clCreateKernel", err);
   
   
   
   return thekernel;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

void getProgramBinaries(cl_program theprogram)
{
   size_t sizes[10];
   char **binaries;
   int i = 0;
   int numdevices=0;
   int totalsize=0;
   int err;

   err = clGetProgramInfo(theprogram, CL_PROGRAM_BINARY_SIZES,
               10*sizeof(size_t), sizes, NULL);
   errchk("clGetProgramInfo: CL_PROGRAM_BINARY_SIZES", err);

   while (sizes[i] > 0)
   {
      numdevices++;
      totalsize += sizes[i];
      printf("size of binary %d is %ld\n", i, sizes[i]);
      i++;
   }

   printf("total num of devices program associated with: %d\n", numdevices);

   binaries = new char *[numdevices];

   for (i=0; i<numdevices; i++)
      binaries[i] = new char[sizes[i]];

   err = clGetProgramInfo(theprogram, CL_PROGRAM_BINARIES,
              numdevices*sizeof(char *), binaries, NULL);
   errchk("clGetProgramInfo: CL_PROGRAM_BINARIES", err);
  
   int j;

   for (i=0; i<numdevices; i++)
   {
//      if (sizes[i] > 0)
//         printf("%s\n", binaries[i]);

      for (j=0; j<sizes[i]; j++)
         printf("%c", binaries[i][j]);

      delete [] binaries[i];
   }

   delete [] binaries;
}

