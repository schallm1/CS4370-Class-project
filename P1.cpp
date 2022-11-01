#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <utility>
#include <vector>
#include <map>

using namespace std;

float maxVal(float a, float b);

int main(int argc, char *argv[])
{
    unordered_multimap<float, string> geneMap;
    unordered_map<string, string> expressionMap;
    vector<string> finalGenes;
    map<string, int> finalCtrl;
    map<string, int> finalCase;

    if (argc !=5)
    {
        printf("Error, five args needed");
        exit(1);
    }
    string rowString;
    string gene;

    int missingControl = 0, missingCases=0;
    float controlSum = 0, caseSum = 0;
    int max = 0;
    int num;
    int place;
    float temp = 0;
    float youden = 0;
    fstream file;

    string name = argv[1];
    char f;

    file.open(argv[1], fstream::in);
        if(!file)
        {
            cout << "File did not open properly. Try again.";
            exit(1);
        }

        int genes = atoi(argv[2]);
        const int cases = atoi(argv[3]);
        const int controlCases = atoi(argv[4]);
        const int matrixSize = controlCases+cases;

        int controlData[controlCases];
        int diseaseData[cases];
        int *matrix = new int[genes*matrixSize];

        if(genes<=0 || cases<=0 || controlCases<=0)
        {
            cout << "The number of genes and individuals must each be greater than one.\n";
            exit(1);
        }
        //skip over header information
        for(int i = 0; i < 1; i++)
        {
            for(int j = 0; j < cases + controlCases + 1; j++)
            {
                file >> rowString;
            }
        }

        for (int x = 0; x<genes; x++)
        {
            caseSum=0;
            controlSum = 0;
            missingCases = 0;
            missingControl = 0;
            temp = 0;
            file >> gene;
            for(int y = 0; y<cases; y++)
            {
                file >> rowString;
                if((rowString[0] < '0' || rowString[0]> '1') && rowString[0] != '-')
                {
                    diseaseData[y] = -2;
                    *(matrix + x * matrixSize + y) = atoi(rowString.c_str());
                }
                else
                {
                    diseaseData[y] = atoi(rowString.c_str());
                    *(matrix + x * matrixSize + y) = atoi(rowString.c_str());
                }
            }
            for(int z = 0; z<controlCases; z++)
            {
                file >> rowString;
                if((rowString[0] < '0' || rowString[0]> '1') && rowString[0] != '-')
                {
                    controlData[z] = -2;
                    *(matrix + x * matrixSize + cases + z) = atoi(rowString.c_str());
                }
                else
                {
                    controlData[z] = atoi(rowString.c_str());
                    *(matrix + x * matrixSize + cases + z) = atoi(rowString.c_str());
                }
            }
            for(int b = 0; b<cases; b++)
            {
                if(diseaseData[b]!=0 && diseaseData[b]!=-1 && diseaseData[b]!=1)
                {
                    missingCases++;
                    continue;
                }
                caseSum+= diseaseData[b];
            }
            if(caseSum <= 0)
            {
                if(caseSum == 0)
                    expressionMap[gene] = "Neutral";
                else
                {
                    caseSum = 0-(caseSum/(cases-missingCases));
                    expressionMap[gene] = "Low";
                }
            }
            else
            {
                caseSum = caseSum/(cases-missingCases);
                expressionMap[gene] = "High";
            }

            for(int a = 0; a<controlCases; a++)
            {
                if(controlData[a]!=0 && controlData[a]!=-1 && controlData[a]!=1)
                {
                    missingControl++;
                    continue;
                }
                controlSum+=controlData[a];
            }
            if(controlSum < 0)
            controlSum = 0-(controlSum/(controlCases-missingControl));
            else
            controlSum = controlSum/(controlCases-missingControl);

            temp = caseSum - controlSum;
            if(temp<0)
            {
                temp = 0-temp;
            }
            geneMap.insert(make_pair(temp, gene));

            youden = maxVal(temp, youden);
            if(youden<0)
            {
                youden = 0-youden;
            }
        }
            int buck = geneMap.bucket(youden);
            pair<unordered_multimap<float, string>::iterator, unordered_multimap<float,string>::iterator> p1 = geneMap.equal_range(youden);
            for(; p1.first != p1.second; p1.first++)
            {
                finalGenes.push_back(p1.first->second);
            }
                unordered_map<string, string>::iterator it;
                map<string, int>::iterator iter;
                cout<< "Expression pattern:\n";
                string disease;
                string ctrl;
                for(int k = 0; k < finalGenes.size(); k++)
                {
                    disease = "";
                    ctrl = "";
                    it = expressionMap.find(finalGenes.at(k));
                    f = it->first[1];
                    place = f;
                    place -= 49;
                    cout << it->first << "\t" << it->second <<endl;
                    if(it->second == "Low")
                    {
                        for(int f = 0; f<cases; f++)
                        {
                            if(*(matrix + place*matrixSize +f)== -1)
                            {
                                disease += "D";
                                num = f+1;
                                disease += to_string(num);
                                if(finalCase.find(disease) == finalCase.end())
                                {
                                    finalCase.emplace(disease, 1);
                                }
                                disease = "";
                                
                            }
                        }
                        for(int n = 0; n<controlCases; n++)
                        {
                            if(*(matrix + place * matrixSize + (n+cases))== 1 || *(matrix + place * matrixSize + (n+cases)) == 0 )
                            {
                                ctrl += "C"; 
                                num = n+1;
                                ctrl += to_string(num);
                                if(finalCtrl.find(ctrl) == finalCtrl.end())
                                {
                                    finalCtrl.emplace(ctrl, 1);
                                }
                                ctrl = "";
                            }
                        }
                    }
                    else if(it->second == "High")
                    {
                        for(int g = 0; g<cases; g++)
                        {
                            if(*(matrix + place * matrixSize + g) == 1)
                            {
                                disease += "D";
                                num = g+1;
                                disease += to_string(num);
                                if(finalCase.find(disease) == finalCase.end())
                                {
                                    finalCase.emplace(disease, 1);
                                }
                                disease = "";
                            }
                        }
                        for(int l = 0; l<controlCases; l++)
                        {
                            if(*(matrix + place * matrixSize + (l+cases)) == -1 || *(matrix + place * matrixSize + (l+cases)) == 0 )
                            {
                                ctrl += "C";
                                num = l+1;
                                ctrl += to_string(num);
                                if(finalCtrl.find(ctrl) == finalCtrl.end())
                                {
                                    finalCtrl.emplace(ctrl, 1);
                                }
                                ctrl = "";
                            }
                        }
                    }
                    else if(it->second == "Neutral")
                    {
                        for(int h = 0; h<cases; h++)
                        {
                            if(*(matrix + place * matrixSize + h) == 0)
                            {
                                disease += ("D");
                                num = h+1;
                                disease += to_string(num);
                                if(finalCase.find(disease) == finalCase.end())
                                {
                                    finalCase.emplace(disease, 1);
                                }
                                disease = "";
                            }
                        }
                        for(int m = 0; m<controlCases; m++)
                        {
                            if(*(matrix + place * matrixSize + (m+cases)) == -1 || *(matrix + place * matrixSize + (m+cases)) == 1 )
                            {
                                ctrl += "C";
                                num = m+1;
                                ctrl += to_string(num);
                                if(finalCtrl.find(ctrl) == finalCtrl.end())
                                {
                                    finalCtrl.emplace(ctrl, 1);
                                }
                                ctrl = "";
                            }
                        }
                    }
                }
                cout << "Cases with pattern:" << endl;
                for(iter = finalCase.begin(); iter!=finalCase.end(); iter++)
                {
                    cout << iter->first << " ";
                }
                cout <<endl << "Controls with pattern: "<<endl;
                for(iter = finalCtrl.begin(); iter!=finalCtrl.end(); iter++)
                {
                    cout << iter->first << " ";
                }
                cout << endl << endl << "j = " << youden;
                delete[] matrix;

}

float maxVal(float a, float b)
{
    if(a>=b)
    {  
       return a; 
    }
    else
    return b;
}