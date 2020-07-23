//
// Created by berry on 2020-04-22.
//

#ifndef SOAPNUKE_PEFILTERTMPOUT_H
#define SOAPNUKE_PEFILTERTMPOUT_H


class peFilterTmpOut {
public:
    bool readType;
    bool filter;
    long long lineNumber;
public:
    peFilterTmpOut(bool readType,bool filter,long long lineNumber){
        this->readType=readType;
        this->filter=filter;
        this->lineNumber=lineNumber;
    }
    string toString(){
        ostringstream concat;
        concat<<readType<<"\t"<<filter<<"\t"<<lineNumber;
        return concat.str();
    }
};


#endif //SOAPNUKE_PEFILTERTMPOUT_H
