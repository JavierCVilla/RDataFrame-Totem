#include "TTree.h"
#include "TFile.h"
#include "ROOT/RDataFrame.hxx"
#include <iostream>

bool CompareTree(const std::string &fileA, const std::string &treeA,
                 const std::string &fileB, const std::string &treeB){

    // Check number of branches
    ROOT::RDataFrame dA(treeA, fileA);
    ROOT::RDataFrame dB(treeB, fileB);

    auto branchesA = dA.GetColumnNames();
    auto branchesB = dB.GetColumnNames();

    if(branchesA != branchesB){
        std::cout << "Trees have different branches" << std::endl;

        std::cout << "Branches in " << fileA << std::endl;
        for(auto b : branchesA)
            std::cout << b << std::endl;

        std::cout << "Branches in " << fileB << std::endl;
        for(auto b : branchesB)
            std::cout << b << std::endl;

        return 1;
    }

    // Read second TTree as a TTree friend
    auto f = new TFile(fileA.c_str());
    auto T = (TTree*)f->Get(treeA.c_str());
    std::ostringstream ft;
    ft << "ft=" << treeB;
    T->AddFriend(ft.str().c_str(), fileB.c_str());

    // Create RDF out of both trees
    ROOT::RDataFrame d(*T);

    auto branches = d.GetColumnNames();

    // Custom filter with an AND of all branches
    std::ostringstream filter;
    filter << "true";
    int offset = branches.size() / 2;

    // Consider only branch names
    // not branch.leaf names
    for (int i = 0; i < offset ; i+=2){
        auto firstElem  = branches[i];
        auto secondElem = branches[i+offset];
        // Go forward to equivalent column
        filter << " && " << firstElem << " == " << secondElem;
    }

    std::cout << "\nFilter: " << std::endl;
    std::cout <<  filter.str() << "\n" << std::endl;
    auto r = d.Define("results", filter.str());

    // TBS: auto all = [](bool a, bool b){return a == b;};
    //auto equalTrees = r.Reduce(all, "results", true);
    if( d.Count().GetValue() == r.Sum("results").GetValue()){
        std::cout << "Both trees are equivalent" << std::endl;
        return 0;
    }else{
        std::cout << "Different trees!" << std::endl;
        return 1;
    }
}

int main(int argc, char* argv[]){
    if(argc < 5){
        std::cout << "Usage: rootcompare fileA treenameA fileB treenameB" << std::endl;
        return 1;
    }
    auto fileA = argv[1];
    auto treeA = argv[2];
    auto fileB = argv[3];
    auto treeB = argv[4];

    return CompareTree(fileA, treeA, fileB, treeB);
}
