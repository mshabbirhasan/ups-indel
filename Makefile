CXX = g++ -std=c++11
INCLUDE = -I include

default: all
all: ups_indel ups_generate_redundant_indel_list ups_compare_uvcf_files ups_filter_uvcf_file

ups_indel:
	$(CXX) $(INCLUDE) src/ups_indel.cpp -lpthread -o ups_indel  

ups_generate_redundant_indel_list:
	$(CXX) $(INCLUDE) src/ups_generate_redundant_indel_list.cpp -o ups_generate_redundant_indel_list 

ups_compare_uvcf_files:
	$(CXX) $(INCLUDE) src/ups_compare_uvcf_files.cpp -o ups_compare_uvcf_files

ups_filter_uvcf_file: 
	javac -d . src/GenerateFilteredUVCFFileAfterRemovingRedundantIndel.java 

clean:
	rm -f ups_indel ups_generate_redundant_indel_list ups_compare_uvcf_files GenerateFilteredUVCFFileAfterRemovingRedundantIndel.class

.phony: clean default
