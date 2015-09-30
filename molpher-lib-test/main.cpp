/* 
 * File:   main.cpp
 * Author: sichom
 *
 * Created on September 23, 2015, 10:09 AM
 */

#include <cstdlib>

#include "testing/testing.hpp"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    std::string storage_dir("Results");
    std::string template_file("../templates/test-template.xml");
    run_path_finder(storage_dir, template_file, 2);
    return 0;
}