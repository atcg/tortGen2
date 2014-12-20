//
//  main.cpp
//  cPWP
//
//  Created by Evan McCartney-Melstad on 12/18/14.
//  Copyright (c) 2014 Evan McCartney-Melstad. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <sstream>

using namespace std;

// Function prototypes:
int coverage(const unsigned char readNums[], int locus, int tortNumber, int totalTorts);
double calcMajorAlleleFreq(const unsigned char readNums[], int locus, int tortNumber, int totalTorts);
long double selfPWP(const unsigned char readNums[], int locus, int tortNumber, int totalTorts, int cov, int locWeighting, double majAlleleFreq);
double interPWP(int tortNum, int comparisonTortNum, int locWeighting, double majAlleleFreqTort, double majAlleleFreqComparisonTort);



int main(int argc, const char * argv[]) {
    if (argc != 3) {
        cout << "You must supply two arguments to this program: 1) the binary readcounts file basename, and 2) and integer showing how many consituent files there are. The filenames should look like <basename>.nmil.binary8bitunsigned, where n is the number of the file." << endl;
        return 1;
    }
    
    int numFiles = atoi(argv[2]);
    string fileBaseName = argv[1];
    string fileExtension = "mil.binary8bitunsigned";
    
    long double pairwisePi[272][272] = {0.0}; // Will hold the pairwise pi estimates
    unsigned long long int weightings[272][272] = {0}; // Will hold the full weightings across all loci
    
    for (int fileNum = 1; fileNum <= numFiles; fileNum++) {
        //string fileName = fileBaseName + std::to_string(fileNum) + fileExtension;
        streampos size; // This will hold the number of bytes in the input file currently being worked on
        
        std::stringstream fileNameBuilder;
        fileNameBuilder << fileBaseName << fileNum << fileExtension;
        string fileName = fileNameBuilder.str();    // unnecessary step, but shown to demonstrate how to obtain the string from stringstream.
        
        ifstream file(fileName.c_str(), ios::in|ios::binary|ios::ate);
        if (file.is_open()) {
            
            size = file.tellg();
            unsigned char* readCounts;
            readCounts = new unsigned char[size];
            
            file.seekg (0, ios::beg);
            file.read((char*)readCounts, size);
            file.close();
            
            cout << "The entire contents of " << fileName << " is in memory." << endl;
            cout << "The total number of sites (given 272 tortoises) in " << fileName << " is " << (int)(size/(272*2)) << endl;
            
            int totalLoci = (int)(size / (272*2));
            
            cout << "totalLoci in " << fileName << " = " << totalLoci << endl;
            
            for (int locus = 0; locus < totalLoci; locus++) {
                unsigned int locusCoverages[272];
                long double majorAlleleFreq[272];
                
                for (int tort = 0; tort < 272; tort++) {
                    
                    locusCoverages[tort] = coverage(readCounts, locus, tort, 272);
                    if(locusCoverages[tort] > 0) {
                        majorAlleleFreq[tort] = calcMajorAlleleFreq(readCounts, locus, tort, 272);
                        if(locusCoverages[tort] > 1) {
                            unsigned int locusWeighting = locusCoverages[tort]*(locusCoverages[tort]-1);
                            weightings[tort][tort] = weightings[tort][tort] + locusWeighting;
                            // Add the pwp estimate for the current locus to the whole
                            pairwisePi[tort][tort] = pairwisePi[tort][tort] + selfPWP(readCounts, locus, tort, 272, locusCoverages[tort], locusWeighting, majorAlleleFreq[tort]);
                        }
                        
                        for (int comparisonTort = 0; comparisonTort <= tort; comparisonTort++) {
                            if (locusCoverages[comparisonTort] > 0) {
                                unsigned int locusWeighting = locusCoverages[tort] * locusCoverages[comparisonTort];
                                weightings[tort][comparisonTort] = weightings[tort][comparisonTort] + locusWeighting;
                                pairwisePi[tort][comparisonTort] = pairwisePi[tort][comparisonTort] + interPWP(tort, comparisonTort, locusWeighting, majorAlleleFreq[tort], majorAlleleFreq[comparisonTort]);
                            }
                        }
                    }
                }
            }
            delete[] readCounts;
        } else {
            cout << "Couldn't read in " << fileName << endl;
            return 1;
        }
        
    }
    
    
    ofstream outFile;
    outFile.open("pwpEstimate.txt");
    if (outFile.is_open()) {
        for (int tort = 0; tort < 272; tort++) {
            for (int comparisonTort = 0; comparisonTort <= tort; comparisonTort++) {
                if (weightings[tort][comparisonTort] > 0) {
                    outFile << tort << ":" << comparisonTort << "\t" << pairwisePi[tort][comparisonTort] / weightings[tort][comparisonTort] << endl;
                    cout << "weighting\t" << weightings[tort][comparisonTort] << endl;
                    cout << "pwp\t" << pairwisePi[tort][comparisonTort] / weightings[tort][comparisonTort] << endl;
                } else {
                    outFile << "NA" << endl;
                }
            }
        }
    } else {
        cout << "Unable to open pwpEstimate.txt for writing" << endl;
        return 1; // Failure
    }
    return 0;
}
        




// Helper Functions:
int coverage(const unsigned char readNums[], int locus, int tortNumber, int totalTorts) {
    /*
     Coverage should be calculated by summing the readcounts for the minor and major alleles
     for a tortoise. The major allele is the first column, and the minor allele is the second column.
     So, for the first tortoise (torts[0]), we want readCounts[0] and readCounts[1] for the first locus.
     The first locus will take up the first totalTorts*2 integers in the readCounts array, so 0-551 in the
     case of 276 tortoises. We want the second locus to therefore be readCounts[552] and readCounts[553].
     */
    int cov = (int)readNums[(locus * totalTorts * 2) + tortNumber*2] + (int)readNums[(locus * totalTorts * 2) + (tortNumber*2) + 1];
    return cov;
}


double calcMajorAlleleFreq(const unsigned char readNums[], int locus, int tortNumber, int totalTorts) {
    double maj = (double)readNums[(locus * totalTorts * 2) + tortNumber*2] / ( (double)readNums[(locus * totalTorts * 2) + (tortNumber*2) + 1] + (double)readNums[(locus * totalTorts * 2) + tortNumber*2] );
    return maj;
}


long double selfPWP(const unsigned char readNums[], int locus, int tortNumber, int totalTorts, int cov, int locWeighting, double majAlleleFreq) {
    //PI[N,N] += WW * ( 2 * P[N] * ( D[N] - $(2*N-1) ) / (D[N]-1) ) ;  # prob of difference
    
    double sfPWP = double(locWeighting) * ( 2.0 * double(majAlleleFreq) * (       (double(cov) - double(readNums[(locus*totalTorts*2)+(tortNumber*2)]))    ) / double(cov-1.0));
    return sfPWP;
}



double interPWP(int tortNum, int comparisonTortNum, int locWeighting, double majAlleleFreqTort, double majAlleleFreqComparisonTort ) {
    //PI[N,M] += WW * ( P[N]*(1-P[M]) + P[M]*(1-P[N]) ) ;
    /* For each locus, the maximum pwp value will be 1 * (locCoverage[tort] * locCoverage[comparisonTort].
     This is on the order of 1 * product total bases sequenced per tort, or 2.5*10^19 if each has 5 billion
     mapped read bases return pwp.
     */
    
    double pwp = double(locWeighting) * ( double(majAlleleFreqTort) * (1.0-double(majAlleleFreqComparisonTort)) + double(majAlleleFreqComparisonTort)*(1.0-double(majAlleleFreqTort)));
    if (pwp < 0.0 || pwp > 625.0) {
        cout << "Warning: Unreasonable pwp estimate for tort #s " << tortNum << " and " << comparisonTortNum << ". Locus PWP = " << pwp << ". Locus weighting = " << locWeighting << ". Tort1 major allele frequency = " << majAlleleFreqTort << ". Tort2 major allele frequency = " << majAlleleFreqComparisonTort << endl;
    }
    return pwp;
}