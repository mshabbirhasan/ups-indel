#ifndef utility
#define utility

#include <string>
#include <cstring>
#include <sstream>
#include <vector>

using namespace std;

vector<string> split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

string diff(string str1, string str2){
	//str1 is always longer than str2

	stringstream result;
    int left = 0, right = str1.length()-1;
    for(int j=0; j<str2.length(); j++){
        if(str1.at(j) == str2.at(j)){
            left++;
        }
        else{
            break;
        }
    }

    if(left == str2.length()){ //CACA CA
        result<<str1.substr(left, str1.length()-left);
    }
    else{
        for(int k = str2.length()-1, l = str1.length()-1; k>0; k--, l--){
            if(str1.at(l) == str2.at(k)){
                right--;
                if(k <= left){
                    break;
                }
            }
            else{
                break;
            }
        }

        result<<str1.substr(left,right-left+1);
    }

	return result.str();
}

#endif
