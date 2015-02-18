// SIGMOD Programming Contest 2015 - simple test driver
//
// This code simply passes a test file to a test program
//---------------------------------------------------------------------------
// This is free and unencumbered software released into the public domain.
//
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
//
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// For more information, please refer to <http://unlicense.org/>
//---------------------------------------------------------------------------
#include <iostream>
#include <map>
#include <chrono>
#include <vector>
#include <cassert>
#include <cstdint>
#include <fcntl.h>
#include <unistd.h>
#include <sys/sendfile.h>
#include <sys/socket.h>
//---------------------------------------------------------------------------
using namespace std;
//---------------------------------------------------------------------------
static int openFile(const char* name)
{
   int fd=open(name,O_RDONLY|O_CLOEXEC);
   if (fd<0) {
      cerr << "unable to open " << name << endl;
      exit(1);
   }
   return fd;
}
//---------------------------------------------------------------------------
static int openClient(const char* program)
{
   int fds[2];
   if (socketpair(AF_UNIX,SOCK_STREAM,0,fds)<0) {
      perror("testdriver");
      exit(1);
   }
   int pid=vfork();
   if (pid<0) {
      perror("testdriver");
      exit(1);
   }
   if (pid==0) {
      close(fds[0]);
      dup2(fds[1],0);
      dup2(fds[1],1);
      close(fds[1]);
      execl("/bin/sh","sh","-c",program,nullptr);
      perror("testdriver");
      _exit(127);
   }
   close(fds[1]);
   return fds[0];
}
//---------------------------------------------------------------------------
static void readFully(int fd,void* target,unsigned len,const char* message)
{
   while (len) {
      int res=read(fd,target,len);
      if (res<1) {
         cout << message << endl;
         exit(1);
      }
      target=static_cast<char*>(target)+res;
      len-=res;
   }
}
//---------------------------------------------------------------------------
int main(int argc,char* argv[])
{
   if (argc!=3) {
      cout << "usage: " << argv[0] << " client testfile" << endl;
      return 1;
   }
   int file=openFile(argv[2]);
   int client=openClient(argv[1]);
   vector<char> buffer,buffer2;
   unsigned clientOffset=0;
   auto start = std::chrono::system_clock::now();
   while (true) {
      uint32_t len;
      int res=read(file,&len,sizeof(uint32_t));
      if (!res) break;
      if (res!=sizeof(uint32_t)) {
         cout << "warning: test file corruption!" << endl;
         return 1;
      }
      if (sendfile(client,file,nullptr,len)!=len) {
         cout << "sending to client failed" << endl;
         return 1;
      }
      if (read(file,&len,sizeof(uint32_t))!=sizeof(uint32_t)) {
         cout << "warning: test file corruption!" << endl;
         return 1;
      }
      if (len>buffer.size()) {
         buffer.resize(len);
         buffer2.resize(len);
      }
      readFully(file,buffer.data(),len,"test file corruption");
      readFully(client,buffer2.data(),len,"reading result from client failed");
      for (uint32_t index=0;index!=len;++index) {
         if (buffer[index]!=buffer2[index]) {
            cout << "wrong validation result for validation " << (index+clientOffset) << ", got " << static_cast<int>(buffer2[index]) << ", expected " << static_cast<int>(buffer[index]) << endl;
            return 1;
         }
      }
      clientOffset+=len;
   }
   auto end = std::chrono::system_clock::now();
   cout << "{{" << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "}}" << endl;
}
//---------------------------------------------------------------------------
