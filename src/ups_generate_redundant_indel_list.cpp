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
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
// BUT NOT LIMITED TO, PROCUREMENT
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
#include <map>
#include <list>
#include<algorithm>
#include "utility.h"

using namespace std;

int main(int argc, char const ** argv){
    if (argc != 3){
        std::cerr << "USAGE: ups_generate_redundant_indel_list UPS_INDEL_VCF_FILE.UVCF OUTPUT_FILE_NAME.\n";
        return 1;
    }

	string line;
    ifstream myInputFile(argv[1]); //input uvcf file
    ofstream myOutputFile(argv[2]); //output text file

    std::map<string, std::list<string> > myHashMap;

    if (myInputFile.is_open()){
        while (getline(myInputFile, line)){
			if(line.at(0)!='#' && line.at(0)!='\n'){
				vector<string> tokens = split(line, '\t');
                if(tokens.at(3).length() != tokens.at(4).length()){
                    string pattern = tokens.at(7); //the whole pattern in the UVCF file Column 7

                    string key = pattern;
                    string value = tokens.at(2);

                    if(myHashMap.count(key) > 0){
                        list<string> valueList = myHashMap.at(key);
                        if(!(std::find(valueList.begin(), valueList.end(),value)!= valueList.end())){
                            myHashMap[key].push_back(value);
                        }
                    }
                    else{
                        myHashMap[key].push_back(value);
                    }
                }
			}
        }
        myInputFile.close();
    }
    else{
        std::cerr << "ERROR: Can't open the input file.\n";
        return 1;
    }

    if(myOutputFile.is_open()){
        for(std::map<string, std::list<string> >::iterator it = myHashMap.begin(); it!= myHashMap.end(); it++){
            std::list<string> myList = it->second;
            if(myList.size() > 1){
                int i = 1;
                myList.sort();
                myOutputFile<<"[";
                for(std::list<string>::iterator it = myList.begin(); it!= myList.end(); it++, i++){
                    myOutputFile<<*it;
                    if(i<myList.size()){
                        myOutputFile<<", ";
                    }
                }
                myOutputFile<<"]"<<endl;
            }
        }
        myOutputFile.close();
    }
    else{
        std::cerr << "ERROR: File output error.\n";
        return 1;
    }

    return 0;
}
