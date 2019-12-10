/*
MIT License

Copyright (c) 2019 Hamada Laboratory, and Ryoga Ishida

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <string>
#include <stdlib.h>
#include <iostream>
#include "RaptRanker.hpp"

using namespace std;

int main(int argc, char* argv[]){
    //DO NOT use printf/scanf
    std::ios::sync_with_stdio(false);

    if (argc !=2) {
        cout << "Error: The number of argument is invalid." << endl;
        return(0);
    }
    
    string parameter_file = argv[1];
    RaptRanker raptranker(parameter_file);
    raptranker.Run();

    
    return(0);
}
