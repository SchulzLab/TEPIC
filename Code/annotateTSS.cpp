//
// Created by Fabian Kern on 27.04.16.
//

#include "annotateTSS.hpp"

gtfdictionary readGTF(char *file) {

    gtfdictionary dict;
    std::ifstream infile(file, std::ifstream::in);
    std::string line;

    while (std::getline(infile, line)) {
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens.size() >= 9) {
            if (tokens[2] == "gene") {
                std::string key = tokens[0];
                replaceAll(key, "chr", "");

                hashitem second;
                std::string first = tokens[9];

                if (tokens[6] == "+") {
                    second = hashitem(key, std::stoi(tokens[3]));
                }
                else {
                    second = hashitem(key, std::stoi(tokens[4]));
                }

                std::pair<std::string, hashitem> entry(first, second);
                dict.insert(entry);
            }
        }

    }

    return dict;
}


ocdictionary readOC_Region(char *file) {

    ocdictionary dict;
    std::ifstream infile(file, std::ifstream::in);
    std::string line;

    while (std::getline(infile, line)) {
        std::vector<std::string> tokens;
        tokens = split(line, ' ');
        tokens = split(tokens[0], ':');

        if (tokens.size() >= 2) {
            std::vector<std::string> subtokens = split(tokens[1], '-');
            std::string key;
            key = tokens[0];
            replaceAll(key, "chr", "");

            std::map<std::string, std::vector<std::pair<int, int> > >::iterator it;
            it = dict.find(key);

            std::pair<int, int> pair(std::stoi(subtokens[0]), std::stoi(subtokens[1]));
            if (it != dict.end()) {
                it->second.push_back(pair);
            }
            else {
                std::vector<std::pair<int, int> > list;
                list.push_back(pair);
                std::pair<std::string, std::vector<std::pair<int, int> > > entry(key, list);
                dict.insert(entry);
            }
        }

    }

    return dict;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

void replaceAll(std::string &str, const std::string &from, const std::string &to) {
    if (from.empty())
        return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
}
