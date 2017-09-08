import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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

public class GenerateFilteredUVCFFileAfterRemovingRedundantIndel {

    public static void main(String args[]) throws FileNotFoundException, IOException {
        if (args.length != 2) {
            System.out.println("USAGE: java GenerateFilteredUVCFFileAfterRemovingRedundantIndel UVCF_FILE REDUNDANT_INDEL_LIST");
            System.exit(1);
        }
        String uvcfFile = args[0];
        String redundantIndelList = args[1];

        FileReader redundantIndelListReader = new FileReader(redundantIndelList);
        BufferedReader redundantIndelListBufferedReader = new BufferedReader(redundantIndelListReader);

        FileReader uvcfFileReader = new FileReader(uvcfFile);
        BufferedReader uvcfBufferedReader = new BufferedReader(uvcfFileReader);

        String outputFileName = uvcfFile.substring(0, uvcfFile.lastIndexOf(".")) + "_filtered.uvcf";
        FileWriter outputFileWriter = new FileWriter(outputFileName);

        try {
            Map<String, List<String>> redundantIndelMap = new HashMap<String, List<String>>(); // if indels r1, r2, and r3 are redundant then r1 is the key and <r2, r3> is the value
            Set<String> processedIndels = new HashSet<String>(); // keeps the indels and their redundant ones so that they are not printed in the output

            String line;
            while ((line = redundantIndelListBufferedReader.readLine()) != null) {
                String lineWithoutSquareBracket1 = line.replace("[", "");
                String lineWithoutSquareBracket2 = lineWithoutSquareBracket1.replace("]", "");

                String lineArray[] = lineWithoutSquareBracket2.split(",");
                List<String> redundantIndels = new ArrayList<String>();
                for (int i = 1; i < lineArray.length; i++) {
                    redundantIndels.add(lineArray[i].trim());
                }
                redundantIndelMap.put(lineArray[0].trim(), redundantIndels);
            }

            while ((line = uvcfBufferedReader.readLine()) != null) {
                if(line.startsWith("#")){
						outputFileWriter.write(line+"\n");
				} else{
					String lineArray[] = line.split("\t");
					String indelId = lineArray[2];

					if (redundantIndelMap.containsKey(indelId)) {
						outputFileWriter.write(line);
						processedIndels.add(indelId);

						List<String> indelList = redundantIndelMap.get(indelId);
						
						outputFileWriter.write(";Redundant Indel Ids: ");
						for (String indel : indelList) {     
							outputFileWriter.write(indel + " ");
							processedIndels.add(indel);
						}
						outputFileWriter.write("\n");
					} else if (!processedIndels.contains(indelId)) {
						outputFileWriter.write(line + "\n");
						processedIndels.add(indelId);
					}
				}
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        } finally {
            outputFileWriter.close();
            redundantIndelListReader.close();
            redundantIndelListBufferedReader.close();
            uvcfFileReader.close();
            uvcfBufferedReader.close();
        }
    }
}
