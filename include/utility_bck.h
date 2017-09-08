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
	stringstream result;
	if(str1.length() > str2.length()){
        if(str1.length()>=2 && str2.length()>=2){
            int lengthDifference = str1.length() - str2.length();

            if(str1.substr(lengthDifference, str1.length()) == str2){ // GTGA GA
                result<<str1.substr(0, lengthDifference);
                return result.str();
            }

            else if(str1.substr(0, str2.length()) == str2){
                result<<str1.substr(str2.length(), str1.length());
                return result.str();
            }

            if(str1.at(0) == str2.at(0) && str2.at(str2.length()-1) == str1.at(str1.length()-1)){
                if(str2.at(0) == str1.at(str1.length()-2)){//TATTG TG
                    result<<str1.substr(1,str1.length()-2);
                    return result.str();
                }
                else{ //AGTC ATC
                    int left = 0, right = str1.length()-1;
                    for(int j=0; j<str2.length(); j++){
                        if(str1.at(j) == str2.at(j)){
                            left++;
                        }
                        else{
                            break;
                        }
                    }
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
                    return result.str();
                }
            }
        }

        int i;
		bool isSpecial = false; // this variable is used to represent the situations where there is any indel in the middle of the REF or ALT seq.
		for(i = 1; i<(int)str2.length(); i++){
			if(str1.at(i)!=str2.at(i) && str1.at(i+1) == str2.at(i)){
                isSpecial = true;
				result<<str1.at(i);
			}
		}

        if(i<(int)str1.length() && isSpecial){
            for(int j = i+1; j<(int)str1.length(); j++){
                result<<str1.at(j);
            }
		}
		else if(i<(int)str1.length() && !isSpecial){
            for(int j = i; j<(int)str1.length(); j++){
                result<<str1.at(j);
            }
		}

		isSpecial = false;
	}
	else if(str2.length() > str1.length()){
        if(str1.length() == 2 && str2.length() == 4){
            if((str1 == str2.substr(0,str1.length())) && (str1 == str2.substr(str1.length(), str1.length()))){ // TA TATA
                result<<str1;
                return result.str();
            }
            else if((str1 == str2.substr(0,str1.length())) && (str1 != str2.substr(str1.length(), str1.length()))){ //TA TATG
                result<<str2.substr(str1.length(), str1.length());
                return result.str();
            }
            else if((str1 != str2.substr(0,str1.length())) && (str1 == str2.substr(str1.length(), str1.length()))){ //TG TATG
                result<<str2.substr(0,str1.length());
                return result.str();
            }
            else if(str1.at(0) == str2.at(0) && str1.at(str1.length()-1) != str2.at(1) && str1.at(str1.length()-1) == str2.at(str2.length()-1)){ //CG CTAG
                result<<str2.substr(1,str2.length()-2);
                return result.str();
            }
        }

        if(str1.compare(str2.substr(1,str2.length()))==0){
            result<<str2.at(0);
            return result.str();
        }

        int lengthDifference = str2.length() - str1.length();
        if((str2.substr(lengthDifference, str2.length()) == str1) && str1.length()>1 ){
            result<<str2.substr(0, lengthDifference);
            return result.str();
        }

        if(str1.length()>2 && str2.length()>2){
            if(str1.at(0) == str2.at(0) && str2.at(str2.length()-1) == str1.at(str1.length()-1)){
                if(str1 != str2.substr(0,str1.length())){
                    if(str1.at(0) == str2.at(str2.length()-2)){
                        result<<str2.substr(1,str2.length()-2);
                        return result.str();
                    }
                    else{ //ATC AGGTC
                        int left = 0, right = str2.length()-1;
                        for(int j=0; j<str1.length(); j++){
                            if(str2.at(j) == str1.at(j)){
                                left++;
                            }
                            else{
                                break;
                            }
                        }

                        for(int k = str1.length()-1, l = str2.length()-1; k>=0; k--, l--){
                            if(str2.at(l) == str1.at(k)){
                                right--;
                                if(k <= left){
                                    break;
                                }
                            }
                            else{
                                break;
                            }
                        }

                        result<<str2.substr(left,right-left+1);
                        return result.str();
                    }
                }
            }
        }

        int i;
        bool isSpecial = false;
        for(i = 1; i<(int)str1.length(); i++){
            if(str2.at(i)!=str1.at(i) && str2.at(i+1) == str1.at(i)){
                isSpecial = true;
                result<<str2.at(i);
            }
        }

        if(i<(int)str2.length() && isSpecial){
            for(int j = i+1; j<(int)str2.length(); j++){
                result<<str2.at(j);
            }
        }

        else if(i<(int)str2.length() && !isSpecial){
            for(int j = i; j<(int)str2.length(); j++){
                result<<str2.at(j);
            }
        }

        isSpecial = false;
	}

	return result.str();
}

#endif
