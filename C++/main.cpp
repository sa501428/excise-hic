/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <iostream>
#include <string>
#include "straw.h"
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 4) {
        cerr << "Incorrect arguments" << endl;
        cerr << "Usage: excise <hicFile> <resolution> <out_short_mnd>" << endl;
        exit(1);
    }

    string filename = argv[1];
    string resolutionString = argv[2];
    int32_t resolution = stoi(resolutionString);

    FILE *fp;
    fp = fopen(argv[3],"w");

    // defaults for extraction
    string norm = "NONE";
    string unit = "BP";
    string matrixType = "observed";

    vector<chromosome> chromosomes = getChromosomes(filename);
    for (int i = 1; i < chromosomes.size(); i++) {
        for (int j = i; j < chromosomes.size(); j++) {
            vector<contactRecord> records = straw(matrixType, norm, filename,
                                                  chromosomes[i].name, chromosomes[j].name, unit, resolution);
            for (contactRecord record : records) {
                auto realCounts = static_cast<int32_t>(record.counts);
                fprintf(fp, "%s %d %s %d %d\n", chromosomes[i].name.c_str(), record.binX,
                        chromosomes[j].name.c_str(), record.binY, realCounts);
            }
        }
    }

    fclose(fp);
    return 0;
}
