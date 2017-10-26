//
//  main.cpp
//  Rosalind
//
//  Created by Mateusz Adamczak on 11/10/2017.
//  Copyright © 2017 Mateusz Adamczak. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <map>


std::map<std::string, std::string> RNATABLE = {
    {"UUU", "F"}, {"CUU", "L"}, {"AUU", "I"}, {"GUU", "V"},
    {"UUC", "F"}, {"CUC", "L"}, {"AUC", "I"}, {"GUC", "V"},
    {"UUA", "L"}, {"CUA", "L"}, {"AUA", "I"}, {"GUA", "V"},
    {"UUG", "L"}, {"CUG", "L"}, {"AUG", "M"}, {"GUG", "V"},
    {"UCU", "S"}, {"CCU", "P"}, {"ACU", "T"}, {"GCU", "A"},
    {"UCC", "S"}, {"CCC", "P"}, {"ACC", "T"}, {"GCC", "A"},
    {"UCA", "S"}, {"CCA", "P"}, {"ACA", "T"}, {"GCA", "A"},
    {"UCG", "S"}, {"CCG", "P"}, {"ACG", "T"}, {"GCG", "A"},
    {"UAU", "Y"}, {"CAU", "H"}, {"AAU", "N"}, {"GAU", "D"},
    {"UAC", "Y"}, {"CAC", "H"}, {"AAC", "N"}, {"GAC", "D"},
    {"UAA", "Stop"}, {"CAA", "Q"}, {"AAA", "K"}, {"GAA", "E"},
    {"UAG", "Stop"}, {"CAG", "Q"}, {"AAG", "K"}, {"GAG", "E"},
    {"UGU", "C"}, {"CGU", "R"}, {"AGU", "S"}, {"GGU", "G"},
    {"UGC", "C"}, {"CGC", "R"}, {"AGC", "S"}, {"GGC", "G"},
    {"UGA", "Stop"}, {"CGA", "R"}, {"AGA", "R"}, {"GGA", "G"},
    {"UGG", "W"}, {"CGG", "R"}, {"AGG", "R"}, {"GGG", "G"}
};

std::string readDataFromFile(const char * filePath){
    std::ifstream fileObject(filePath);
    std::stringstream buffer;

    buffer << fileObject.rdbuf();
    
    return buffer.str();
}

std::vector<std::string> readFileByLines(const char * filePath){
    std::vector<std::string> lines;
    std::ifstream input(filePath);
    
    for (std::string line; getline(input, line, '\n');){
        lines.push_back(line);
    }
    
    return lines;
}



std::vector<std::tuple<std::string, std::string>> readGCContentFile(const char * filePath){
    std::ifstream input(filePath);
    std::string identifier;
    std::string content;
    
    std::vector<std::tuple<std::string, std::string>> gcContentVector;
    
    for (std::string line; getline(input, line);){
        if (line.substr(0, 10) == ">Rosalind_"){
            if (!content.empty()){
                gcContentVector.push_back(std::make_tuple(identifier, content));
                content = "";
            }
            identifier = line;
        }
        else{
            content += line;
        }
    }
    
    if (!content.empty()){
        gcContentVector.push_back(std::make_tuple(identifier, content));
    }
    

    return gcContentVector;
}

int countNucleotide(std::string dna, const char * nucleotide){
    int count = 0;
    for (int i = 0; i < dna.length(); i ++){
        if(dna.at(i) == *nucleotide){
            count ++;
        }
    }
    return count;
}

std::string transcribeToRna(std::string dna){
    const char * thymine = "T";
    
    for (int i = 0; i < dna.length(); i++){
        const char current = dna.at(i);
        if (current == *thymine){
            dna.replace(i, 1 ,"U");
        }
    }
    return dna;
}


std::string complementaryDnaStrand(std::string dna){
    const char * ad = "A";
    const char * ct = "C";
    const char * gu = "G";
    const char * th = "T";
    //write complementary
    std::string complementary;
    
    for(int i =0; i < dna.length(); i++){
        if (dna.at(i) == *ad){
            complementary.append("T");
        }
        else if(dna.at(i) == *ct){
            complementary.append("G");

        }
        else if(dna.at(i) == *gu){
            complementary.append("C");
            
        }
        else if(dna.at(i) == *th){
            complementary.append("A");
            
        }
    }
    std::reverse(complementary.begin(), complementary.end());
    
    //reverse it
    return complementary;
}


float countGCContent(std::string dnaFragment){
    int gcount = countNucleotide(dnaFragment, "G");
    int ccount = countNucleotide(dnaFragment, "C");
    int acount = countNucleotide(dnaFragment, "A");
    int tcount = countNucleotide(dnaFragment, "T");
    //std::cout << gcount << ", " << ccount << ", " << acount << ", "<< tcount << "\n";
    
    return (gcount + ccount) / float(gcount + ccount + acount + tcount);
}

void findHighestGCContent(std::vector<std::tuple<std::string, std::string>> dnaData){
    std::string name;
    float content = 0.0f;
    
    for(std::tuple<std::string, std::string> tpl : dnaData){
        std::string tempName = std::get<0>(tpl);
        float tempContent = countGCContent(std::get<1>(tpl));
        if (tempContent > content){
            name = tempName;
            content = tempContent;
        }
    }
    std::cout << name << "\n" << content * 100 << std::endl;

}

int countPointMutations(std::string first, std::string second){
    if (first.length() != second.length()){
        return int(first.length());
    }
    
    int count = 0;
    
    for (size_t i = 0; i < first.length(); i++){
        if (first.at(i) != second.at(i)){
            count++;
        }
    }
    
    return count;
}

std::string translateRnaIntoProtein(std::string rna){
    std::string protein;
    
    for(size_t i = 0; i < rna.length(); i=i+3){
        protein += RNATABLE[rna.substr(i, 3)];
        }
    std::cout << protein << std::endl;

    return protein;
}


void findMotifInDna(std::string dna, std::string motif){
    size_t index = dna.find(motif);
    while (index < dna.length()){
        std::cout << index + 1 << " ";
        index = dna.find(motif, index + 1);
    }
    std::cout << std::endl;
}

void FindMostCommon(std::vector<std::tuple<std::string, std::string>> dnaData){
    std::string consensus;
    size_t length = std::get<1>(dnaData.at(0)).length();
    int ac[length];
    int gc[length];
    int cc[length];
    int tc[length];
    
    for (size_t i = 0; i < length; i++){
        std::string current;
        for (std::tuple<std::string, std::string> dna : dnaData){
            current += std::get<1>(dna).at(i);
        }
        
        const char * temp = "A";
        int maxhit = 0;
        
        ac[i] = countNucleotide(current, "A");
        maxhit = ac[i];
        
        gc[i] = countNucleotide(current, "G");
        
        if (gc[i] > maxhit){
            temp = "G";
            maxhit = gc[i];
        }
        
        cc[i] = countNucleotide(current, "C");
        
        if (cc[i] > maxhit){
            temp = "C";
            maxhit = cc[i];
        }
        tc[i] = countNucleotide(current, "T");
        if (tc[i] > maxhit){
            temp = "T";
            maxhit = tc[i];
        }
        
        consensus += temp;
 
    }
    
    std::cout << consensus << std::endl;
    
    std::cout << "A: ";
    for (size_t j = 0; j < length; j++){
        std::cout << ac[j] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "C: ";
    for (size_t j = 0; j < length; j++){
        std::cout << cc[j] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "G: ";
    for (size_t j = 0; j < length; j++){
        std::cout << gc[j] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "T: ";
    for (size_t j = 0; j < length; j++){
        std::cout << tc[j] << " ";
    }
    std::cout << std::endl;
    
}


std::string findSharedMotif(std::vector<std::tuple<std::string, std::string>> dnaData){
    std::string longest;
    int hits = 0;
    std::string firstOne = std::get<1>(dnaData.at(0));
    
    for (size_t i = 0; i < firstOne.length(); i++){
        for (size_t j = firstOne.length(); j > i ; j--){
            std::string lookFor = firstOne.substr(i, j - i);
            for (std::tuple<std::string, std::string> tpl : dnaData){
                if (std::get<1>(tpl).find(lookFor) < std::get<1>(tpl).length()){
                    hits++;
                }
            }
            
            if (hits == dnaData.size() && lookFor.length() > longest.length() ){
                longest = lookFor;
            }
            else{
                hits = 0;
            }
        }
    }
    return longest;
}



int main(int argc, const char * argv[]) {
    /*
     Finding a Shared Motif
     */
    
    std::vector<std::tuple<std::string, std::string>> dnaData = readGCContentFile("/Users/mateusz.adamczak/Downloads/rosalind_lcsm.txt");
    std::cout << findSharedMotif(dnaData) << std::endl;
    
    
    /*
     Consensus and Profile
     
    std::vector<std::tuple<std::string, std::string>> dnaData = readGCContentFile("/Users/mateusz.adamczak/Downloads/rosalind_cons.txt");
    
    FindMostCommon(dnaData);
     
     */
    
    
    
    /*
     Finding a Motif in DNA
     std::string dna = readDataFromFile("/Users/mateusz.adamczak/Downloads/rosalind_subs-2.txt");
    
    findMotifInDna(dna, "CTTATTACT");
     */
    
    /*
     Translating RNA into Protein
    std::string rna = readDataFromFile("/Users/mateusz.adamczak/Downloads/rosalind_prot-3.txt");
    
    translateRnaIntoProtein(rna);
     
     */
    
    
    /*
     Counting Point Mutations
    std::vector<std::string> data = readFileByLines("/Users/mateusz.adamczak/Downloads/rosalind_hamm.txt");
    std::cout << countPointMutations(data.at(0), data.at(1)) << std::endl;
    */
    
    
    /*
     Computing GC Content
     
    std::vector<std::tuple<std::string, std::string>> dnaData = readGCContentFile("/Users/mateusz.adamczak/Downloads/rosalind_gc-4.txt");
    findHighestGCContent(dnaData);
     */

    
    
    /*
     Complementing a Strand of DNA
     
    std::string dna = readDataFromFile("/Users/mateusz.adamczak/Downloads/rosalind_revc.txt");
    std::cout << complementaryDnaStrand(dna) << "\n";
     */
    
    
    /*
     Transcribing DNA into RNA
     
    std::string dna = readDataFromFile("/Users/mateusz.adamczak/Downloads/rosalind_rna.txt");
    std::cout << transcribeToRna(dna);
     */

    /*
     Counting DNA Nucleotides
    std::string dna = readDataFromFile("/Users/mateusz.adamczak/Downloads/rosalind_dna.txt");
    std::cout << countNucleotide(dna, "A")  << " ";
    std::cout << countNucleotide(dna, "C")  << " ";
    std::cout << countNucleotide(dna, "G")  << " ";
    std::cout << countNucleotide(dna, "T")  << "\n";
     */
    
    
    
    return 0;
}
