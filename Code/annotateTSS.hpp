//
// Created by Fabian Kern on 27.04.16.
//

#ifndef TEPIC_HIC_ANNOTATETSS_HPP
#define TEPIC_HIC_ANNOTATETSS_HPP


#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

typedef std::pair <std::string, int> hashitem;
typedef std::map <std::string, hashitem> gtfdictionary;
typedef std::map <std::string, std::vector<std::pair<int, int> > > ocdictionary;

gtfdictionary readGTF(char* file);
std::vector<std::string> split(const std::string &s, char delim);
void replaceAll(std::string& str, const std::string& from, const std::string& to);

#endif //TEPIC_HIC_ANNOTATETSS_HPP
