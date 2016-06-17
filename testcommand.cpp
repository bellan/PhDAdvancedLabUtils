#include <iostream>
#include <stdio.h>
#include <stdexcept>
#include <string>

using namespace std;

std::string exec(const char* cmd) {
    char buffer[128];
    std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try {
        while (!feof(pipe)) {
            if (fgets(buffer, 128, pipe) != NULL)
                result += buffer;
        }
    } catch (...) {
        pclose(pipe);
        throw;
    }
    pclose(pipe);
    return result;
}

int main() {
  
  cout << exec("ls") << endl;
  exec("paplay /usr/share/sounds/ubuntu/stereo/bell.ogg");
}

