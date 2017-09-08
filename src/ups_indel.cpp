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
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
// WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Mohammad Shabbir Hasan (shabbir5@vt.edu)
// ==========================================================================

#include <cstdlib>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <fcntl.h>
#include <cstdio>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

#include <pthread.h>

#include "utility.h"

#include <iostream>
#include <fstream>

using namespace std;

bool isHorizontalDecompositionEnabled = true;
string reference_file_name;
string input_file_name;
string input_file_base_name;
string output_file_name;

inline bool fileExists (const string& fileName) {
    if (FILE *file = fopen(fileName.c_str(), "r")) {
        fclose(file);
        return true;
    }
    return false;
}

bool isComplexVariant(string r, string a){
    if(r.length() > 1 && a.length() > 1){
        if(r.at(1) != a.at(1)){
            return true;
        }
    }
    return false;
}

string getLeftCircularPermutation(string * str) {
    string s = str -> substr(1, str -> length());
    string s2(1, str -> at(0));

    s += s2;

    return s;
}

string getRightCircularPermutation(string * str) {
    string s(1, str -> at(str -> length() - 1));
    string s2 = str -> substr(0, str -> length() - 1);

    s += s2;

    return s;
}

void generateEir(string chrom, string pos, string id, string refer, string alt,
    string qual, string filter, string info, string * sequence,
    unsigned long int position, string pattern, string type,
    FILE * fout) {
    if (type == "I") {
        string x;
        string x_prime;
        string str1;
        string str2;
        unsigned long int ir, il;

        int length_difference = alt.length() - refer.length();

        if (refer.length() > 1 && refer == alt.substr(0, refer.length())) {
            position += refer.length() - 1;
        } else if (refer.length() > 2 &&
            refer == alt.substr(length_difference, alt.length())) {
            position -= 1;
        }

        /*********Extend to the right *************/
        x = pattern;
        x_prime = getLeftCircularPermutation( & x);
        ir = position - 1; // because in VCF file position starts from 1 but in

        if ((ir + 1) < sequence -> length()) {
            string r(1, sequence -> at(ir + 1));

            str1 = x + r;
            str2 = r + x_prime;

            while (str1.compare(str2) == 0) {
                x = x_prime;
                ir++;
                if ((ir + 1) >= sequence -> length()) {
                    ir--;
                    break;
                }
                r = sequence -> at(ir + 1);
                x_prime = getLeftCircularPermutation( & x);

                str1 = x + r;
                str2 = r + x_prime;
            }
        } else {
            ir--;
        }

        /*********Extend to the left *************/
        x = pattern;
        x_prime = getRightCircularPermutation( &x);
        il = position - 1; // because in VCF file position starts from 1 but seqan
        // starts string index from 0.
        if (il > 0) {
            string l(1, sequence -> at(il));

            str1 = l + x;
            str2 = x_prime + l;

            while (str1.compare(str2) == 0) {
                x = x_prime;
                il--;
                if (il < 0) {
                    break;
                }
                l = sequence -> at(il);
                x_prime = getRightCircularPermutation( &x);

                str1 = l + x;
                str2 = x_prime + l;
            }
        }

        fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s%lu%s%lu%s\t%s\n",
            chrom.c_str(), pos.c_str(), id.c_str(), refer.c_str(), alt.c_str(),
            qual.c_str(), filter.c_str(), "+", x.c_str(), "[", (il + 2), " - ",
            (ir + 2), "]", info.c_str());
    } else if (type == "D") {
        string x;
        string x_prime;
        string str1;
        string str2;
        unsigned long int ir, il;

        int deletionLength = refer.length() - alt.length();

        if ((alt.length() > 1 && alt == refer.substr(0, alt.length()))) {
            if (alt != refer.substr(deletionLength, refer.length())) { // TCA TC
                position += alt.length() - 1;
            } else {
                if (alt.length() <= 2) { // AAAAAAAAA AA
                    position += alt.length() - 1;
                } else { // AAACAAA AAA
                    position += deletionLength - 1;
                }
            }
        } else if (alt.length() == 1 && alt == refer.substr(0, alt.length()) &&
            alt == refer.substr(refer.length() - 1)) { // CAC C
            if (deletionLength == 1) { // CC C
                position -= 1;
            } else { // CAC C
                position -= alt.length() + 1;
            }
        } else if (alt == refer.substr(deletionLength, refer.length()) &&
            refer.length() != 4 && alt.length() != 2 &&
            alt != refer.substr(0, deletionLength - 1)) { // TATAC TAC
            position -= 1;
        }

        /*********Extend to the right *************/

        x = pattern;

        int length = strlen(pattern.c_str());
        x_prime = getLeftCircularPermutation( & x);

        ir = position + 1;
        if (ir < sequence -> length()) {
            string substring = sequence -> substr(ir, length);

            while (substring.compare(x_prime) == 0) {
                ir++;
                x = x_prime;
                x_prime = getLeftCircularPermutation( & x);
                substring = sequence -> substr(ir, length);
            }
        }

        /*********Extend to the left *************/
        x = pattern;
        x_prime = getRightCircularPermutation( & x);
        il = position - 1;
        if (il >= 0) {
            string l(1, sequence -> at(il));

            str1 = l + x;
            str2 = x_prime + l;

            while (str1.compare(str2) == 0) {
                x = x_prime;
                il--;
                if (il < 0) {
                    break;
                }
                l = sequence -> at(il);
                x_prime = getRightCircularPermutation( & x);

                str1 = l + x;
                str2 = x_prime + l;
            }
        }

        fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s%s%lu%s%lu%s\t%s\n",
            chrom.c_str(), pos.c_str(), id.c_str(), refer.c_str(), alt.c_str(),
            qual.c_str(), filter.c_str(), "-", x.c_str(), "[", (il + 2), " - ",
            ir, "]", info.c_str());
    } else if (type == "N/A") {
        fprintf(fout, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", chrom.c_str(),
            pos.c_str(), id.c_str(), refer.c_str(), alt.c_str(), qual.c_str(),
            filter.c_str(), "N/A[]", info.c_str());
    }
}

int genenateUPSCoordinate(string referenceFileName, string inputFileName,
    string outputFileName) {

    ifstream fin(inputFileName);
    FILE * fout = fopen((outputFileName + ".uvcf").c_str(), "w");

    string type;

    bool isFirstTime = true;
    string chromosomeNumber;
    seqan::FaiIndex faiIndex;
    unsigned idx = 0;
    seqan::Dna5String sequenceInfix;
    stringstream seqStr;
    string reference;

    int ups_rs_id = 1;

    if (fout != NULL && fin.is_open()) {
        string line;
        while (getline(fin, line)) {
            if (line.at(0) == '#') {
                if (line.at(1) == '#') {
                    fprintf(fout, "%s\n", line.c_str()); // only the header of the VCF file
                } else {
                    fprintf(
                        fout, "%s\n",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tUPS-COORDINATE\tINFO");
                }
            } else if (line.at(0) != '#' && line.at(0) != '\n') {
                vector < string > tokens = split(line, '\t');
                std::size_t found;
                string pattern;
                string chrom = tokens.at(0); // CHROM

                if (isFirstTime) {
                    chromosomeNumber = chrom;

                    // Try to load index and create on the fly if necessary.
                    if (seqan::read(faiIndex, referenceFileName.c_str()) != 0) {
                        if (build(faiIndex, referenceFileName.c_str()) != 0) {
                            std::cerr << "ERROR: Index could not be loaded or built.\n";
                            return 1;
                        }

                        // Name is stored from when reading.
                        if (write(faiIndex) != 0) {
                            std::cerr << "ERROR: Index could not be written do disk.\n";
                            return 1;
                        }
                    }

                    // Translate sequence name to index.
                    if (!getIdByName(faiIndex, chromosomeNumber, idx)) {
                        std::cerr << "ERROR: Index does not know about sequence " << referenceFileName << ".\n";
                        return 1;
                    }

                    // Finally, get infix of sequence.
                    if (readSequence(sequenceInfix, faiIndex, idx) != 0) {
                        std::cerr << "ERROR: Could not load infix.\n";
                        return 1;
                    }

                    seqStr << sequenceInfix;
                    reference = seqStr.str();

                    isFirstTime = false;
                } else if (chrom != chromosomeNumber) {
                    string chromosomeNumber = chrom;

                    // Translate sequence name to index.
                    unsigned idx = 0;
                    if (!getIdByName(faiIndex, chromosomeNumber, idx)) {
                        std::cerr << "ERROR: Index does not know about sequence " << referenceFileName << ".\n";
                        return 1;
                    }

                    // Finally, get infix of sequence.
                    if (readSequence(sequenceInfix, faiIndex, idx) != 0) {
                        std::cerr << "ERROR: Could not load infix.\n";
                        return 1;
                    }

                    stringstream newSeqStr;
                    newSeqStr << sequenceInfix;

                    reference = newSeqStr.str();
                }

                string pos = tokens.at(1); // POS
                string id = tokens.at(2); // ID

                // in some cases, there is no id in the ID column of vcf, in that case
                // assign an id
                if (id == ".") {
                    id = std::to_string(ups_rs_id);
                    ups_rs_id++;
                }
                string refer = tokens.at(3); // REF
                string alt = " "; // ALT will be filled later
                string qual = tokens.at(5); // QUAL
                string filter = tokens.at(6); // FILTER
                string info = tokens.at(7); // INFO

                found = tokens.at(4).find(','); // some rows have more than one indel in the alt column
                if (found != string::npos) {
                    vector < string > alternative_indel_token =
                        split(tokens.at(4).c_str(), ',');
                    for (int i = 0; i < alternative_indel_token.size(); i++) {
                        if ((alternative_indel_token.at(i)).at(0) != '<') {
                            alt = alternative_indel_token.at(i);
                            if (refer.length() < alt.length()) {
                                type = "I"; // insertion
                                if(isHorizontalDecompositionEnabled){
                                    pattern = diff(alt,refer); // pattern means difference between ref sequence and alt sequence in the VCF file
                                }
                                else{
                                    if(isComplexVariant(refer, alt)){
                                        type = "N/A";
                                        pattern = "[]";
                                    }
                                    else{
                                        pattern = diff(alt, refer);
                                    }
                                }
                            } else if (refer.length() > alt.length()) {
                                type = "D"; // deletion
                                if(isHorizontalDecompositionEnabled){
                                    pattern = diff(refer,alt); // pattern means difference between ref sequence and alt sequence in the VCF file
                                }
                                else{
                                    if(isComplexVariant(refer, alt)){
                                        type = "N/A";
                                        pattern = "[]";
                                    }
                                    else{
                                        pattern = diff(refer, alt);
                                    }
                                }
                            } else if (refer.length() == alt.length()) {
                                type = "N/A"; // when REF and ALT have same length, such as SNPs and MNPs
                                pattern = "[]";
                            }

                            generateEir(chrom, pos, id, refer, alt, qual, filter, info, & reference, strtoul(pos.c_str(), NULL, 0),
                                pattern, type, fout);
                        }
                    }
                } else {
                    alt = tokens.at(4); // ALT

                    // indel without any comma
                    if (refer.length() < alt.length()) {
                        type = "I";
                        if(isHorizontalDecompositionEnabled){
                            pattern = diff(alt, refer);
                        }
                        else{
                            if(isComplexVariant(refer, alt)){
                                type = "N/A";
                                pattern = "[]";
                            }
                            else{
                                pattern = diff(alt, refer);
                            }
                        }
                    } else if (refer.length() > alt.length()) {
                        type = "D";
                        if(isHorizontalDecompositionEnabled){
                            pattern = diff(refer, alt);
                        }
                        else{
                            if(isComplexVariant(refer, alt)){
                                type = "N/A";
                                pattern = "[]";
                            }
                            else{
                                pattern = diff(refer, alt);
                            }
                        }
                    }

                    // when REF and ALT have same length, such as SNPs and MNPs
                    else if (refer.length() == alt.length()) {
                        type = "N/A";
                        pattern = "[]";
                    }

                    generateEir(chrom, pos, id, refer, alt, qual, filter, info, & reference, strtoul(pos.c_str(), NULL, 0),
                        pattern, type, fout);
                }
            }
        }
        fin.close();
        fclose(fout);
    } else {
        std::cerr << "ERROR: File operation error.\n";
        return 1;
    }

    return 0;
}

int splitInput(string inputFileName) {
    string command = "java -jar ext/SnpSift.jar split " + inputFileName;
    return system(command.c_str());
}

void * executeUPSIndel(void * threadIndex) {
    long chrNum;
    chrNum = (long) threadIndex;

    int flag = -999;
    string chr = "";

    if (chrNum >= 1 && chrNum <= 22) {
        chr = to_string(chrNum);
    } else if (chrNum == 23) {
        chr = "X";
    } else if (chrNum == 24) {
        chr = "Y";
    } else if (chrNum == 25) {
        chr = "MT";
    }

    string vcfFileName = input_file_base_name + "." + chr + ".vcf";

    if(fileExists(vcfFileName)){
        flag = genenateUPSCoordinate(reference_file_name,
        vcfFileName, output_file_name + "." + chr);

        if (flag) {
            exit(-1);
        }
    }
}

int merge_outputs() {
    bool hasHeaderPrinted = false;
    ofstream fout(output_file_name + ".uvcf");
    if (fout.is_open()) {
        for (int chrNum = 1; chrNum <= 25; chrNum++) {
            string chr = "";
            if (chrNum >= 1 && chrNum <= 22) {
                chr = to_string(chrNum);
            } else if (chrNum == 23) {
                chr = "X";
            } else if (chrNum == 24) {
                chr = "Y";
            } else if (chrNum == 25) {
                chr = "MT";
            }

            string chromosome_wise_uvcf_file_name =
                output_file_name + "." + chr + ".uvcf";
            string line;

            ifstream uvcfFile(chromosome_wise_uvcf_file_name);
            if (fileExists(chromosome_wise_uvcf_file_name)) {
                if(uvcfFile.is_open()){
                    while (getline(uvcfFile, line)) {
                        if (line.at(0) == '#' && !hasHeaderPrinted) {
                            fout << line << endl;
                        } else if (line.at(0) != '#') {
                            fout << line << endl;
                        }
                    }
                    uvcfFile.close();
                } else {
                    std::cerr << "ERROR: Could not read from UVCF file: " << chromosome_wise_uvcf_file_name << ".\n";
                    return 1;
                }

                hasHeaderPrinted = true;
            }
        }

        fout.close();
    } else {
        std::cerr << "ERROR: Could not write output: " << output_file_name << ".uvcf.\n";
        return 1;
    }

    return 0;
}

int cleanup() {
    string remove_temp_vcf_files_command =
        "rm " + input_file_base_name + ".*.vcf";

	if (system(remove_temp_vcf_files_command.c_str()) < 0)
		return 1;

    string remove_temp_uvcf_files_command =
        "rm " + output_file_name + ".*.uvcf";

	if (system(remove_temp_uvcf_files_command.c_str()) < 0)
		return 1;

    return 0;
}

int main(int argc, char const* *argv) {
    if (argc != 5) {
        std::cerr << "USAGE: ups_indel REFERENCE_FILE.fa VCF_FILE.vcf OUTPUT_FILE_NAME -hd=true.\n";
        return 1;
    }

    reference_file_name = argv[1];
    input_file_name = argv[2];
    output_file_name = argv[3];
    string horizontal_decomposition_flag = argv[4];

    if(horizontal_decomposition_flag == "-hd=false"){
        isHorizontalDecompositionEnabled = false;
    }

    int splitFlag = splitInput(input_file_name);
    if (splitFlag < 0) {
      return 1;
    }

    pthread_t threads[25]; // 24 chromosomes and MT

    string input_file_base_name_with_extension =
        input_file_name.substr(0, input_file_name.find_last_of(".vcf") + 1);
    string input_file_name_str(input_file_base_name_with_extension);

    input_file_base_name =
        input_file_name_str.substr(0, input_file_name_str.find_last_of("."));

    for (int k = 0; k < 25; k++) {
        int thread_creation_flag = pthread_create( & threads[k], NULL, executeUPSIndel, (void * )(intptr_t)(k + 1));

        if (thread_creation_flag) {
            std::cerr << "ERROR: Unable to create thread, " << thread_creation_flag << ".\n";
            exit(-1);
        }
    }

    for (int l = 0; l < 25; l++) {
        pthread_join(threads[l], NULL);
    }

    if (merge_outputs()) {
        return 1;
    }

    if (cleanup()) {
      return 1;
    }

    return 0;
}
