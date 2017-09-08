// ==========================================================================
//                                 UPS-indel
// ==========================================================================
// Copyright (c) 2015, Virginia Tech
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL MOHAMMAD SHABBIR HASAN, XIAOWEI WU, 
// LAYNE T. WATSON, LIQING ZHANG OR VIRGINIA TECH BE LIABLE FOR ANY DIRECT, 
// INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Mohammad Shabbir Hasan (shabbir5@vt.edu)
// ==========================================================================

#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>
#include "utility.h"

using namespace std;

int main(int argc, char const ** argv){
    if (argc != 4){
        std::cerr << "USAGE: ups_compare_uvcf_files FIRST_UVCF_FILE.UVCF SECOND_UVCF_FILE.UVCF OUTPUT_FILE_NAME\n";
        return 1;
    }

	string firstUVCF_line, secondUVCF_line;

    ifstream firstUVCFFile(argv[1]); //first uvcf file
    ifstream secondUVCFFile(argv[2]); //second uvcf file
    ofstream outputFile(argv[3]); //output file

    vector<string> first;
    vector<string> second;

    if (firstUVCFFile.is_open()){
        while (getline(firstUVCFFile, firstUVCF_line)){
			if(firstUVCF_line.at(0)!='#' && firstUVCF_line.at(0)!='\n'){
				vector<string> tokens = split(firstUVCF_line, '\t');
                string value = tokens.at(7); // UPS-Coordinate column

                first.push_back(value);
			}
        }
        firstUVCFFile.close();
    }
    else{
        std::cerr << "ERROR: Can't open the first input file: "<<argv[1]<<".\n";
        return 1;
    }

    if(secondUVCFFile.is_open()){
        while (getline(secondUVCFFile, secondUVCF_line)){
			if(secondUVCF_line.at(0)!='#' && secondUVCF_line.at(0)!='\n'){
				vector<string> tokens = split(secondUVCF_line, '\t');
                string value = tokens.at(7); // UPS-Coordinate column

                second.push_back(value);
			}
        }
        secondUVCFFile.close();
    }
    else{
        std::cerr << "ERROR: Can't open the second input file: "<<argv[2]<<".\n";
        return 1;
    }

    vector<string> intersection;
    vector<string> in_first_not_in_second;
    vector<string> in_second_not_in_first;

    sort(first.begin(), first.end());
    sort(second.begin(), second.end());

    set_intersection(first.begin(), first.end(), second.begin(), second.end(), back_inserter(intersection));
    set_difference(first.begin(), first.end(), second.begin(), second.end(), back_inserter(in_first_not_in_second));
    set_difference(second.begin(), second.end(), first.begin(), first.end(), back_inserter(in_second_not_in_first));


    if(outputFile.is_open()){
        outputFile<<"Number of Common Indels in "<<argv[1]<<" and "<<argv[2]<<" : "<<intersection.size()<<endl;
        outputFile<<"Number of Indels in "<<argv[1]<<" but not in "<<argv[2]<<" : "<<in_first_not_in_second.size()<<endl;
        outputFile<<"Number of Indels in "<<argv[2]<<" but not in "<<argv[1]<<" : "<<in_second_not_in_first.size()<<endl;

        outputFile.close();
    }
    else{
        std::cerr << "ERROR: Can't open the output file: "<<argv[3]<<".\n";
        return 1;
    }

    return 0;
}
