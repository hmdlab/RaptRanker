//
//  main.cpp
//  
//
//  Created by adinnovaiton on 2017/02/10.
//
//

#include <string>
#include <stdlib.h>
#include <iostream>
#include "RaptRanker.hpp"

using namespace std;

int main(int argc, char* argv[]){
    //入出力高速化のおまじない
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
