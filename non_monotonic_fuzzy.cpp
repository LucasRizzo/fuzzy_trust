#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath> 
#include <assert.h>
#include <string>
#include <string.h>
#include <iterator>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <limits>
#include "math.h"

using namespace std;

#define EPSILON 0.0001

struct Point2D
{
    double x;
    double y;
};

struct fuzzy_t {
    Point2D p1;
    Point2D p2;
    Point2D p3;
    double mean;
    double sd;
    double mdegree;
    double crispActive;
    double weight;
    string shape;
    string name;
    string level;
    // Flag to indicate whether the degree has already been reduced by attackers
    bool reduced = false;
};

struct rule_t {
    vector<fuzzy_t> fuzzy_variables;
    fuzzy_t fuzzy_index;
    double mdegree;
    double rebuttalMDegree;
    string label;
};

struct attack_t {
    string source;
    string target;
    bool rebuttal;
};

struct mitigating_t {
    vector<string> source;
    vector<string> sourceLevel;
    string target;
    bool composed;
};


// Graphs intersections
double lineY (Point2D p1, Point2D p2, double x);
double normalY(double mean, double sd, double x);
double normalXLeft(double mean, double sd, double y);
double normalXRight(double mean, double sd, double y, double precision);
double lineX (Point2D p1, Point2D p2, double y);
double phi(double x);
double intersectionY(Point2D p1, Point2D p2, Point2D p3, Point2D p4);
double intersectionX(Point2D p1, Point2D p2, Point2D p3, Point2D p4);

// Set methods
void setDegree(double input, fuzzy_t& fuzzy_variable);
void setFuzzyIndex(string level, vector<fuzzy_t>& index, fuzzy_t& fuzzy_index);
void setFuzzyVariable(string name, string level, vector<fuzzy_t>& fuzzy, fuzzy_t& fuzzy_variable);
void addFuzzyToRule(string name, string level, vector<fuzzy_t>& fuzzy, rule_t& rule);
void setRuleDegree(rule_t& rule);

// Print methods
void printFuzzy(vector<fuzzy_t>&);
void printRules(vector<rule_t>& rules, vector<attack_t>& attacks, vector<mitigating_t>& mitigatings);
void printData(vector<vector<string>>& data);
void printPolygon(vector<Point2D> polygon);
void printRulesDegree(vector<rule_t>& rules);

// Read methods
void readFuzzy(vector<fuzzy_t>&, string file);
void readRules(vector<fuzzy_t>& fuzzy_variables, vector<fuzzy_t>& index, vector<rule_t>& rules, vector<attack_t>& attacks, vector<mitigating_t>& mitigatings, string file);
void readData(vector<vector<string>>& data, string file);

// Fuzzy methods
void fuzzyInference(vector<rule_t>& rules, vector<attack_t>& attacks, vector<mitigating_t>& mitigatings, vector<fuzzy_t>& fuzzy_variables, vector<fuzzy_t>& index, vector<string>& data, vector<string>& header, string outPutModel);
Point2D compute2DPolygonCentroid(vector<Point2D>& vertices);
vector<Point2D> definePolygonBiomarkers(vector<double>& maxIndexDegree, vector<fuzzy_t>& index);
vector<Point2D> definePolygonTrust_Linear(vector<double>& maxIndexDegree, vector<fuzzy_t>& index);
vector<Point2D> definePolygonMWL(vector<double>& maxIndexDegree, vector<fuzzy_t>& index);
vector<Point2D> definePolygonMWL_New(vector<double>& maxIndexDegree, vector<fuzzy_t>& index);
vector<Point2D> definePolygonMWL_HCI(vector<double>& maxIndexDegree, vector<fuzzy_t>& index);
vector<Point2D> definePolygonMWL_Normal(vector<double>& maxIndexDegree, vector<fuzzy_t>& index);
vector<Point2D> definePolygonTrust_Normal(vector<double>& maxIndexDegree, vector<fuzzy_t>& index);
void solveAttacks(attack_t& current, vector<rule_t>& rules, vector<attack_t>& attacks, vector<fuzzy_t>& fuzzy_variables);
void solveMitigatings(vector<rule_t>& rules, mitigating_t& mitigating, vector<fuzzy_t>& fuzzy_variables);
void solveRebuttals(vector<rule_t>& rules, attack_t& current);
void findLeaves(attack_t& current, int currentIndex, vector<attack_t>& attacks, vector<rule_t>& rules, vector<int>& leaves);
double meanOfMax (vector<Point2D> inferenceGraph);

// Weight methods
void applyWeights(vector<rule_t>& rules);
void normalizeWeights(vector<fuzzy_t>& fuzzy_variables);

// Get methods
double getValue(vector<string>& data, vector<string>& header, string name);
void getTargetSourceRuleIndex(int& targetIndex, int& sourceIndex, attack_t& attack, vector<rule_t>& rules);
double aggregateRules(vector<rule_t>& rules, string level);
double getWeight(vector<string>& data, vector<string>& header, string name);

// Remove methods
void removeVariable(vector<rule_t>& rules, vector<attack_t>& attacks, string name);

// String operations
int stringOcurrences(string s, string target);
int number_of_lines(string file);

bool inVector(vector<int>& vector, int size, int elem);

// Compile
// g++ -o a -std=c++11 fuzzy.cpp -g


// ------ PARAMETERS ------
// Set fuzzy logic
// Options: Zadeh, Godels, Lukasiewicz, Product
//      Zadeh: AND = min, OR = max
//      Godels: AND = min, OR = max
//      Lukasiewicz: AND = max(a + b - 1, 0), OR = min(a + b, 1)
//      Product: AND = a.b OR = a + b - ab
string fuzzy_logic = "Zadeh";
// Options: MWL, BIO, Trust
string application = "Trust";

bool DEBUG = false;

// Fuzzification graph
bool PolygonBio = false;
bool PolygonMWL = false;
bool PolygonMWL_NEW = false;
bool PolygonMWL_HCI = false;
bool PolygonMWL_Normal = false;
bool PolygonTrust_Linear = true;
bool PolygonTrust_Gaussian = false;

bool WEIGHTS = false;

// ------------------------

int main(int argc, char *argv[]) {

    // Define fuzzy variables
    vector<fuzzy_t> fuzzy_variables;
    vector<fuzzy_t> index;
    vector<vector<string>> data;
    vector<vector<string>> models;
    vector<rule_t> rules;
    vector<attack_t> attacks;
    vector<mitigating_t> mitigatings;

    int truePolygons = 0;

    if (PolygonBio) {
        truePolygons++;
    }

    if (PolygonMWL) {
        truePolygons++;
    }

    if (PolygonMWL_NEW) {
        truePolygons++;
    }

    if (PolygonMWL_HCI) {
        truePolygons++;
    }

    if (PolygonMWL_Normal) {
        truePolygons++;
    }

    if (PolygonTrust_Linear) {
        truePolygons++;
    }

    if (PolygonTrust_Gaussian) {
        truePolygons++;
    }

    if (truePolygons != 1) {
        cout << "Error! Choose one fuzzification polygon!";
        return 0;
    }

    //printFuzzy(fuzzy_variables);
    //printFuzzy(index);
    //printData(data);
    //printRules(rules, attacks, mitigatings);

    if (application == "Trust") {
        
        models.resize(number_of_lines("trust/models.txt"));

        ifstream infile("trust/models.txt");
        string line;

        string delimiter = ",";
        size_t pos = 0;
        string token;

        int models_i = 0;

        while (getline(infile, line)) {
            while ((pos = line.find(delimiter)) != string::npos) {
                token = line.substr(0, pos);
                models[models_i].push_back(token);
                line.erase(0, pos + delimiter.length());
            }

            models[models_i].push_back(line);
            models_i++;
        }

        for (int m = 0; m < models.size(); m++) {

            fuzzy_logic = models[m][1];

            if (models[m][2] == "No") {
                WEIGHTS = false;
            } else {
                WEIGHTS = true;
            }

            if (models[m][3] == "Linear") {
                //readFuzzy(fuzzy_variables, "trust/linear/Kb1 02mar19/parameters.txt");
                //readFuzzy(index, "trust/linear/Kb1 02mar19/trust_index.txt");
                
                readFuzzy(fuzzy_variables, "trust/linear/Kb3 03sep19/parameters.txt");
                readFuzzy(index, "trust/linear/Kb3 03sep19/trust_index.txt");
                
                //readData(data, "trust/italian.csv");
                readData(data, "trust/portuguese.csv");
                
                readRules(fuzzy_variables, index, rules, attacks, mitigatings, "trust/linear/Kb3 03sep19/rules.txt");
                //readRules(fuzzy_variables, index, rules, attacks, mitigatings, "trust/linear/Kb1 02mar19/rules.txt");
                
                //printRules(rules, attacks, mitigatings);

                PolygonTrust_Linear = true;
                PolygonTrust_Gaussian = false;

            } else {
                readFuzzy(fuzzy_variables, "trust/gaussian/Kb3 03sep19/parameters.txt");
                readFuzzy(index, "trust/gaussian/Kb3 03sep19/trust_index.txt");
                
                //readFuzzy(fuzzy_variables, "trust/gaussian/Kb1 02mar19/parameters.txt");
                //readFuzzy(index, "trust/gaussian/Kb1 02mar19/trust_index.txt");
                
                //readData(data, "trust/italian.csv");
                readData(data, "trust/portuguese.csv");
                
                readRules(fuzzy_variables, index, rules, attacks, mitigatings, "trust/gaussian/Kb3 03sep19/rules.txt");
                //readRules(fuzzy_variables, index, rules, attacks, mitigatings, "trust/gaussian/Kb1 02mar19/rules.txt");

                PolygonTrust_Linear = false;
                PolygonTrust_Gaussian = true;
            }

            if (rules.size() == 0) {
                cout << "No rule could be created for associated indexs and variables!";
                return 0;
            }

            cout << "Model: " << models[m][0] << ", " << models[m][1] << ", " << models[m][2] << ", " << models[m][3] << endl;

            // Apply fuzzy inference for each record of data
            // Header is saved at data[0]
            for (int i = 1; i < data.size(); i++) {
                
            //for (int i = 8290; i < 8291; i++) {
            //    cout << "ID: " << i << endl;
                fuzzyInference(rules, attacks, mitigatings, fuzzy_variables, index, data[i], data[0], models[m][0]);
            }

            fuzzy_variables.clear();
            index.clear();
            data.clear();
            rules.clear();
            attacks.clear();
            mitigatings.clear();
        }
    }

    // *** Gambiarra
    /*vector<vector<string>> data2;
    readData(data2, "fuzzy_dataset_100_female.csv");
    vector<rule_t> rules2;
    vector<attack_t> attacks2;
    readRules(fuzzy_variables, index, rules2, attacks2, "bio_rules_5.txt");
    if (rules2.size() == 0) {
        cout << "No rule could be created for associated indexs and variables!";
        return 0;
    }
    for (int i = 1; i < data2.size(); i++) {
        fuzzyInference(rules2, attacks2, fuzzy_variables, index, data2[i], data2[0]);
    }*/

    return 0;
}


void printPolygon(vector<Point2D> polygon) {

    for (int i = 0; i < polygon.size(); i++) {
        cout << polygon[i].x << " " << polygon[i].y << endl;
    }
}

void setRuleDegree(rule_t& rule) {

    double mdegree;

    bool zadehOrGodels = fuzzy_logic.compare("Zadeh") == 0 || fuzzy_logic.compare("Godels") == 0;
    bool lukasiewicz = fuzzy_logic.compare("Lukasiewicz") == 0;
    bool product = fuzzy_logic.compare("Product") == 0;

    if (zadehOrGodels) {
        mdegree = 1;
    }

    if (lukasiewicz) {
        mdegree = 0;
    }

    for (int fuzzy_i = 0; fuzzy_i < rule.fuzzy_variables.size(); fuzzy_i++) {
        if (zadehOrGodels) {
            // Taking the min cause it is assumed that 
            // premissed are connected only by AND
            if (rule.fuzzy_variables[fuzzy_i].mdegree < mdegree) {
                mdegree = rule.fuzzy_variables[fuzzy_i].mdegree;
            }
        } else if (lukasiewicz) {
            mdegree += rule.fuzzy_variables[fuzzy_i].mdegree;
        } else if (product) {
            if (fuzzy_i == 0) {
                mdegree = rule.fuzzy_variables[fuzzy_i].mdegree;
            } else {
                mdegree *= rule.fuzzy_variables[fuzzy_i].mdegree;
            }
        }
    }

    if (lukasiewicz) {
        mdegree = mdegree - (rule.fuzzy_variables.size() - 1);

        if (mdegree < 0) {
            mdegree = 0;
        }
    }

    if (DEBUG) {
        cout << rule.label << ": degree: " << mdegree << " - size: " << rule.fuzzy_variables.size() << endl;
    }
    rule.mdegree = mdegree;
}

void fuzzyInference(vector<rule_t>& rules,
                    vector<attack_t>& attacks,
                    vector<mitigating_t>& mitigatings,
                    vector<fuzzy_t>& fuzzy_variables,
                    vector<fuzzy_t>& index,
                    vector<string>& data,
                    vector<string>& header,
                    string outPutModel) {

    for (int i = 0; i < fuzzy_variables.size(); i++) {

        double value = getValue(data, header, fuzzy_variables[i].name);
        fuzzy_variables[i].weight = getWeight(data, header, fuzzy_variables[i].name);

        if (DEBUG) {
            cout << fuzzy_variables[i].name << ": " << fuzzy_variables[i].level << " = " << value;
        }

        fuzzy_variables[i].mdegree = 0;
        if (value != -999) {
            setDegree(value, fuzzy_variables[i]);
        }

        if (DEBUG) {
            cout << ". Degree: " << fuzzy_variables[i].mdegree << endl;
        }
    }

    normalizeWeights(fuzzy_variables);

    // Set variables degree in the variables inside rules
    for (int rule_i = 0; rule_i < rules.size(); rule_i++) {
        for (int fuzzy_i = 0; fuzzy_i < rules[rule_i].fuzzy_variables.size(); fuzzy_i++) {
            for (int variable_i = 0; variable_i < fuzzy_variables.size(); variable_i++) {
                if (rules[rule_i].fuzzy_variables[fuzzy_i].name == fuzzy_variables[variable_i].name && 
                    rules[rule_i].fuzzy_variables[fuzzy_i].level == fuzzy_variables[variable_i].level) {

                    rules[rule_i].fuzzy_variables[fuzzy_i].mdegree = fuzzy_variables[variable_i].mdegree;
                    rules[rule_i].fuzzy_variables[fuzzy_i].weight = fuzzy_variables[variable_i].weight;
                    break;
                }

                if (variable_i == fuzzy_variables.size() - 1) {
                    cout << "Error! Rule antecedent with no degree!" << endl;
                }
            }
        }
    }

    // Set rule antescedent degree before attacks
    for (int i = 0; i < rules.size(); i++) {
        setRuleDegree(rules[i]);
    }
    
    //cout << "Set rule degree" << endl;
    //printRulesDegree(rules);

    // Solve mitigating attacks
    for (int i = 0; i < mitigatings.size(); i++) {
        solveMitigatings(rules, mitigatings[i], fuzzy_variables);
    }
    
    //cout << "Solve mitigatings" << endl;
    //printRulesDegree(rules);

    vector<int> leaves;
    for (int i = 0; i < attacks.size(); i++) {
        findLeaves(attacks[i], i, attacks, rules, leaves);
    }

    // Solve attacks between rules
    for (int i = 0; i < leaves.size(); i++) {
        solveAttacks(attacks[leaves[i]], rules, attacks, fuzzy_variables);
    }
    
    //cout << "Solve attacks" << endl;
    //printRulesDegree(rules);

    for (int rule_i = 0; rule_i < rules.size(); rule_i++) {
        // After solving non rebuttal attacks, the mdegree of the rules
        // to solve rebuttals will be the mdegree set so far.
        // Rebuttals are applied simultaneously, so it is necessary
        // to have the original value when solving more than one rebuttal.
        rules[rule_i].rebuttalMDegree = rules[rule_i].mdegree;
    }

    for (int i = 0; i < attacks.size(); i++) {
        solveRebuttals(rules, attacks[i]);
    }

    if (DEBUG) {
        cout << "Degrees before weights" << endl;
        printRulesDegree(rules);
    }

    // After solving all the conflicts rule weights will be applied by multiplying
    // the final degree by the rule weight. Since what is being imported is
    // feature weights, the rule weight will be the maximum weight of its features.
    if (WEIGHTS) {
        applyWeights(rules);
    }

    if (DEBUG) {
        cout << "Degrees after weights" << endl;
        printRulesDegree(rules);
    }

    vector<double> maxIndexDegree (index.size());

    bool greaterThanZero = false;
    double max = 0;
    int indexMax = 0;
    for (int i = 0; i < index.size(); i++) {
        maxIndexDegree[i] = aggregateRules(rules, index[i].level);
        if (DEBUG) {
            cout << "max " << i << ": " << maxIndexDegree[i] << endl;
        }

        if (maxIndexDegree[i] > EPSILON) {
            greaterThanZero = true;
        }

        if (max < maxIndexDegree[i]) {
            max = maxIndexDegree[i];
            indexMax = i;
        }
    }
    
    std::ofstream outfile;
    outfile.open(outPutModel, std::ios_base::app);

    if (! greaterThanZero) {
        //cout << "Centroid is null\n";
        outfile << "null, null\n";

    } else {

        vector<Point2D> inferenceGraph;

        if(PolygonBio) {
            inferenceGraph = definePolygonBiomarkers(maxIndexDegree, index);
        }

        if (PolygonMWL) {
            inferenceGraph = definePolygonMWL(maxIndexDegree, index);
        }

        if (PolygonMWL_NEW) {
            inferenceGraph = definePolygonMWL_New(maxIndexDegree, index);
        }

        if (PolygonMWL_HCI) {
            inferenceGraph = definePolygonMWL_HCI(maxIndexDegree, index);
        }

        if (PolygonMWL_Normal) {
            inferenceGraph = definePolygonMWL_Normal(maxIndexDegree, index);
        }

        if (PolygonTrust_Linear) {
            inferenceGraph = definePolygonTrust_Linear(maxIndexDegree, index);
        }

        if (PolygonTrust_Gaussian) {
            inferenceGraph = definePolygonTrust_Normal(maxIndexDegree, index);
        }

        if (DEBUG) {
            printPolygon(inferenceGraph);
        }

        Point2D centroid = compute2DPolygonCentroid(inferenceGraph);

        //cout << data[0] << ", " << centroid.x << "\n";// << ", " << centroid.y << "\n";
        //cout << data[0] << ", " <<  meanOfMax (inferenceGraph) << "\n";
        //cout << data[0] << ", " << centroid.x << ", " << meanOfMax (inferenceGraph) << "\n";
        //cout << centroid.x * 5 << ", " << meanOfMax (inferenceGraph) * 5 << "\n";
        //cout << centroid.x << ", " << centroid.y << ", " << meanOfMax (inferenceGraph) << "\n";

        outfile << centroid.x << ", " << meanOfMax (inferenceGraph) << "\n";

        //cout << centroid.x << ", " << meanOfMax (inferenceGraph) << "\n";
    }
}

double meanOfMax (vector<Point2D> inferenceGraph) {

    // Find max
    double max = inferenceGraph[0].y;
    double index = 0;
    for (int i = 1; i < inferenceGraph.size(); i++) {
        if (max < inferenceGraph[i].y) {
            max = inferenceGraph[i].y;
            index = i;
        }
    }

    // Find ties 
    vector<int> ties;
    for (int i = 0; i < inferenceGraph.size(); i++) {
        if ( abs(inferenceGraph[i].y - inferenceGraph[index].y) < EPSILON && i != index) {
            ties.push_back(i);
        }
    }
    ties.push_back(index);

    double sum = 0;
    for (int i = 0; i < ties.size(); i++) {
        sum += inferenceGraph[ties[i]].x;
    }

    return sum / ties.size();
}



/*
+ Underload                                          Overload
|
XXXXX                     XXXXXXX                     XXXXXXXX
|   XXXX                 X     XXX                  XX
|      XX      Fitting  X        XXX   Fitting     XX
|        XX    Minus   X           X   Plus       XX
|         X            X            X            XX
|         XX          X              X          XX
|          XX         X              XX         X
|           X         X               XX       XX
|            X       X                 XX      X
|            XX      X                   XX    XX
|             XX    X                    XX  XX
|              XX XX                      XXXX
|              XXXX                         XX
|           XXXX  XX                    XXX XXX
+---------XXX------XXX---------------XXXX-----XXXX----------+
0                            50                            100
*/

vector<Point2D> definePolygonMWL_Normal(vector<double>& maxIndexDegree, vector<fuzzy_t>& index) {

    vector<Point2D> polygon;
    vector<Point2D> polygonUnderload;
    vector<Point2D> polygonFittingMinus;
    vector<Point2D> polygonFittingPlus;
    vector<Point2D> polygonOverload;
    int n = maxIndexDegree.size() - 1;

    double PRECISION = 0.01;

    // Underload
    if (maxIndexDegree[0] > EPSILON) {

        polygonUnderload.push_back(Point2D {0, 0});

        if (maxIndexDegree[0] < 1 - EPSILON) {
            // Add line from (0, maxIndexDegree) until intersection with
            // x coordinate of normal curve
            polygonUnderload.push_back(Point2D {0, maxIndexDegree[0]});
            double x = normalXRight(index[0].mean, index[0].sd, maxIndexDegree[0], PRECISION);

            // Add all points of line. This helps to solve intersections
            for (double i = PRECISION; i < x; i += PRECISION) {
                polygonUnderload.push_back(Point2D {i, maxIndexDegree[0]});
            }

            // Add the rest of the normal curve
            for (double i = x; i < index[0].mean + 4*index[0].sd; i += PRECISION) {

                double y = normalY(index[0].mean, index[0].sd, i); 
                if (y < EPSILON) {
                    break;
                }

                polygonUnderload.push_back(Point2D {i, y});
            }
        } else {
            // Degree 1. Add whole normal curve with increments of PRECISION
            double x = index[0].mean;
            do {
                polygonUnderload.push_back(Point2D {x, normalY(index[0].mean, index[0].sd, x)});
                x += PRECISION;
            } while (normalY(index[0].mean, index[0].sd, x) > PRECISION);
        }
    }

    // FittingMinus
    if (maxIndexDegree[1] > EPSILON) {

        double top;

        if (maxIndexDegree[1] < 1 - EPSILON) {
            top = normalXLeft(index[1].mean, index[1].sd, maxIndexDegree[1]);
        } else {
            top = index[1].mean;
        }

        for (double x = index[1].mean - 4*index[1].sd; x <= index[1].mean; x += PRECISION) {
            if (normalY(index[1].mean, index[1].sd, x) < EPSILON) {
                continue;
            }

            if (abs(x - top) < EPSILON) {
                break;
            }

            polygonFittingMinus.push_back(Point2D {x, normalY(index[1].mean, index[1].sd, x)});
        }

        // Make line if maxDegree is smalled than 1
        if (maxIndexDegree[1] < 1 - EPSILON && polygonFittingMinus.size() > 0) {
            // Add all points in line. This helps to solve intersections
            for (double i = polygonFittingMinus[polygonFittingMinus.size() - 1].x + PRECISION; i < index[1].mean; i += PRECISION) {
                polygonFittingMinus.push_back(Point2D {i, maxIndexDegree[1]});
            }
        }

        if (maxIndexDegree[2] < EPSILON) {
            polygonFittingMinus.push_back(Point2D {index[1].mean, 0});
        }
    }

    // FittingPlus
    if (maxIndexDegree[2] > EPSILON) {

        if (maxIndexDegree[1] < EPSILON) {
            polygonFittingPlus.push_back(Point2D {index[2].mean, 0});
        }

        if (maxIndexDegree[2] < 1 - EPSILON) {

            double x = normalXRight(index[2].mean, index[2].sd, maxIndexDegree[2], PRECISION);

            // Draw line
            polygonFittingPlus.push_back(Point2D {index[2].mean, maxIndexDegree[2]});
            // Add all points in line. This helps to solve intersections
            for (double i = index[2].mean + PRECISION; i < x; i += PRECISION) { 
                polygonFittingPlus.push_back(Point2D {i, maxIndexDegree[2]});
            }

            // Draw curve from line
            for (double i = x; i < index[2].mean + 4*index[2].sd; i += PRECISION) {
                if (normalY(index[2].mean, index[2].sd, i) < EPSILON) {
                    break;
                }
                polygonFittingPlus.push_back(Point2D {i, normalY(index[2].mean, index[2].sd, i)});
            }
        } else {
            double x = index[2].mean;

            for (double i = x; i < index[2].mean + 4*index[2].sd; i += PRECISION) {
                if (normalY(index[2].mean, index[2].sd, i) < EPSILON) {
                    break;
                }
                polygonFittingPlus.push_back(Point2D {i, normalY(index[2].mean, index[2].sd, i)});
            }
        }
    }

    // Overload
    if (maxIndexDegree[3] > EPSILON) {

        double top = normalXLeft(index[3].mean, index[3].sd, maxIndexDegree[3]);
        for (double x = index[3].mean - 4*index[3].sd; x < index[3].mean; x += PRECISION) {

            if (x > index[3].mean - EPSILON) {
                break;
            }

            if (normalY(index[3].mean, index[3].sd, x) < EPSILON) {
                continue;
            }

            if (abs(x - top) < EPSILON) {
                break;
            }

            polygonOverload.push_back(Point2D {x, normalY(index[3].mean, index[3].sd, x)});
        }

        // Make line if maxDegree is smalled than 1
        if (maxIndexDegree[3] < 1 - EPSILON && polygonOverload.size() > 0) {
            // Add all points in line. This helps to solve intersections
            for (double i = polygonOverload[polygonOverload.size()-1].x + PRECISION; i < index[3].mean; i += PRECISION) {
                polygonOverload.push_back(Point2D {i, maxIndexDegree[3]});
            }
        }

        polygonOverload.push_back(Point2D {index[3].mean, 0});
    }

    // Remove intersections
    for (int fm = 0; fm < polygonFittingMinus.size(); fm++) {
        for (int un = 0; un < polygonUnderload.size(); un++) {
            if (polygonFittingMinus[fm].y < polygonUnderload[un].y &&
                abs(polygonFittingMinus[fm].x - polygonUnderload[un].x) <= PRECISION) {
                polygonFittingMinus.erase(polygonFittingMinus.begin() + fm);
                fm--;
                break;
            }
        }
    }

    for (int un = 0; un < polygonUnderload.size(); un++) {
        for (int fm = 0; fm < polygonFittingMinus.size(); fm++) {
            if (polygonUnderload[un].y < polygonFittingMinus[fm].y &&
                abs(polygonUnderload[un].x - polygonFittingMinus[fm].x) <= PRECISION) {
                polygonUnderload.erase(polygonUnderload.begin() + un);
                un--;
                break;
            }
        }
    }

    for (int fp = 0; fp < polygonFittingPlus.size(); fp++) {
        for (int ov = 0; ov < polygonOverload.size(); ov++) {
            if (polygonFittingPlus[fp].y < polygonOverload[ov].y && 
                abs(polygonFittingPlus[fp].x - polygonOverload[ov].x) <= PRECISION) {
                polygonFittingPlus.erase(polygonFittingPlus.begin() + fp);
                fp--;
                break;
            }
        }
    }

    for (int ov = 0; ov < polygonOverload.size(); ov++) {
        for (int fp = 0; fp < polygonFittingPlus.size(); fp++) {
            if (polygonOverload[ov].y < polygonFittingPlus[fp].y &&
                abs(polygonOverload[ov].x - polygonFittingPlus[fp].x) <= PRECISION) {
                polygonOverload.erase(polygonOverload.begin() + ov);
                ov--;
                break;
            }
        }
    }

    for (int un = 0; un < polygonUnderload.size(); un++) {
        polygon.push_back(polygonUnderload[un]);
    }

    for (int fm = 0; fm < polygonFittingMinus.size(); fm++) {
        polygon.push_back(polygonFittingMinus[fm]);
    }

    for (int fp = 0; fp < polygonFittingPlus.size(); fp++) {
        polygon.push_back(polygonFittingPlus[fp]);
    }

    for (int ov = 0; ov < polygonOverload.size(); ov++) {
        polygon.push_back(polygonOverload[ov]);
    }

    return polygon;
}

/*
 +Low               Medium low           Medium high          High
 XXXX               XXXXXX               XXXXXXXXX            XX
 X  XXXX           XX    XXX           XX        XX          XX
 |     XX         XX       X           X          X        XX
 |      XX       XX        XX         XX           X      XX
 |       XX      X          XX       XX            XX    XX
 |        XX    X            XX      X              XX   X
 |         X   XX             XX    XX               X   X
 |         XX  X               XX  XX                XX X
 |          XXXX                XXXX                  XXX
 |           XXX                XXXX                 XXX
 |         XXX XXXX          XXXX  XXXX            XXX XX
++-----XXXXX------XXXXXX--XXXX--------XXXX------XXX-----XXXXXXX
0           0.25                0.50                0.75      1
*/

vector<Point2D> definePolygonTrust_Normal(vector<double>& maxIndexDegree, vector<fuzzy_t>& index) {

    vector<Point2D> polygon;
    vector<Point2D> polygonLow (1000);
    vector<Point2D> polygonMediumLow (1000);
    vector<Point2D> polygonMediumHigh (1000);
    vector<Point2D> polygonHigh (1000);
    int n = maxIndexDegree.size() - 1;
    double PRECISION = 0.001;

    // Low
    int pLow = 0;
    if (maxIndexDegree[0] > EPSILON) {

        polygonLow[pLow] = Point2D {0, 0};
        pLow++;

        if (maxIndexDegree[0] < 1 - EPSILON) {
            // Add line from (0, maxIndexDegree) until intersection with
            // x coordinate of normal curve
            polygonLow[pLow] = Point2D {0, maxIndexDegree[0]};
            pLow++;
            double x = normalXRight(index[0].mean, index[0].sd, maxIndexDegree[0], PRECISION);

            // Add all points of line. This helps to solve intersections
            for (double i = PRECISION; i < x; i += PRECISION) {
                polygonLow[pLow] = Point2D {i, maxIndexDegree[0]};
                pLow++;
            }

            // Add the rest of the normal curve
            for (double i = x; i < index[0].mean + 4*index[0].sd; i += PRECISION) {

                double y = normalY(index[0].mean, index[0].sd, i); 
                if (y < EPSILON) {
                    break;
                }

                polygonLow[pLow] = Point2D {i, y};
                pLow++;
            }
        } else {
            // Degree 1. Add whole normal curve with increments of PRECISION
            double x = index[0].mean;
            do {
                polygonLow[pLow] = Point2D {x, normalY(index[0].mean, index[0].sd, x)};
                pLow++;
                x += PRECISION;
            } while (normalY(index[0].mean, index[0].sd, x) > PRECISION);
        }
    }

    // Medium low
    int pMediumLow = 0;
    if (maxIndexDegree[1] > EPSILON) {
        for (double x = index[1].mean - 4*index[1].sd; x <= index[1].mean + 4*index[1].sd; x += PRECISION) {

            double y = normalY(index[1].mean, index[1].sd, x);
            if (y < EPSILON) {
                continue;
            }

            // Line if degree < 1
            if(y >= maxIndexDegree[1]) {
                y = maxIndexDegree[1];
            }

            polygonMediumLow[pMediumLow] = Point2D {x, y};
            pMediumLow++;
        }
    }

    // Medium high
    int pMediumHigh = 0;
    if (maxIndexDegree[2] > EPSILON) {
        for (double x = index[2].mean - 4*index[2].sd; x <= index[2].mean + 4*index[2].sd; x += PRECISION) {

            double y = normalY(index[2].mean, index[2].sd, x);
            if (y < EPSILON) {
                continue;
            }

            // Line if degree < 1
            if(y >= maxIndexDegree[2]) {
                y = maxIndexDegree[2];
            }

            polygonMediumHigh[pMediumHigh] = Point2D {x, y};
            pMediumHigh++;
        }
    }

    // High
    int pHigh = 0;
    if (maxIndexDegree[3] > EPSILON) {

        double top = normalXLeft(index[3].mean, index[3].sd, maxIndexDegree[3]);
        //cout << "TOP: " << top << endl;

        for (double x = index[3].mean - 4*index[3].sd; x < index[3].mean; x += PRECISION) {

            if (x > index[3].mean - EPSILON) {
                break;
            }

            if (normalY(index[3].mean, index[3].sd, x) < EPSILON) {
                continue;
            }

            if (x >= top) {
                break;
            }

            polygonHigh[pHigh] = Point2D {x, normalY(index[3].mean, index[3].sd, x)};
            pHigh++;
        }

        // Make line if maxDegree is smalled than 1
        if (maxIndexDegree[3] < 1 - EPSILON && pHigh > 0) {
            // Add all points in line. This helps to solve intersections
            for (double i = polygonHigh[pHigh-1].x + PRECISION; i < index[3].mean; i += PRECISION) {
                polygonHigh[pHigh] = Point2D {i, maxIndexDegree[3]};
                pHigh++;
            }
        }

        polygonHigh[pHigh] = Point2D {index[3].mean, 0};
        pHigh++;
    }

    // Remove intersections
    vector<int> polygonLowErase (1000);
    int iPLE = 0;
    vector<int> polygonMediumLowErase (1000);
    int iPMLE = 0;
    vector<int> polygonMediumHighErased (1000);
    int iPMHE = 0;
    vector<int> polygonHighErase (1000);
    int iPHE = 0;

    int max = pLow;
    if (pMediumLow > max) {
        max = pMediumLow;
    }
    if (pMediumHigh > max) {
        max = pMediumHigh;
    }
    if (pHigh > max) {
        max = pHigh;
    }

    for (int i = 0; i < max; i++) {
        for (int j = 0; j < max; j++) {
            if (i < pMediumLow && j < pLow &&
                polygonMediumLow[i].y < polygonLow[j].y &&
                abs(polygonMediumLow[i].x - polygonLow[j].x) <= PRECISION) {
                polygonMediumLowErase[iPMLE] = i;
                iPMLE++;
            }

            if (i < pLow && j < pMediumLow && 
                polygonLow[i].y < polygonMediumLow[j].y &&
                abs(polygonLow[i].x - polygonMediumLow[j].x) <= PRECISION) {
                polygonLowErase[iPLE] = i;
                iPLE++;
            }

            if (i < pMediumLow && j < pMediumHigh &&
                polygonMediumLow[i].y < polygonMediumHigh[j].y &&
                abs(polygonMediumLow[i].x - polygonMediumHigh[j].x) <= PRECISION) {
                if (! inVector(polygonMediumLowErase, iPMLE, i)) {
                    polygonMediumLowErase[iPMLE] = i;
                    iPMLE++;
                }
            }

            if (i < pMediumHigh && j < pMediumLow &&
                polygonMediumHigh[i].y < polygonMediumLow[j].y &&
                abs(polygonMediumHigh[i].x - polygonMediumLow[j].x) <= PRECISION) {
                polygonMediumHighErased[iPMHE] = i;
                iPMHE++;
            }

            if (i < pMediumHigh && j < pHigh &&
                polygonMediumHigh[i].y < polygonHigh[j].y && 
                abs(polygonMediumHigh[i].x - polygonHigh[j].x) <= PRECISION) {
                if (! inVector(polygonMediumHighErased, iPMHE, i)) {
                    polygonMediumHighErased[iPMHE] = i;
                    iPMHE++;
                }
            }

            if (i < pHigh && j < pMediumHigh &&
                polygonHigh[i].y < polygonMediumHigh[j].y &&
                abs(polygonHigh[i].x - polygonMediumHigh[j].x) <= PRECISION) {
                polygonHighErase[iPHE] = i;
                iPHE++;
            }
        }
    }

    for (int low = 0; low < pLow; low++) {
        if (! inVector(polygonLowErase, iPLE, low)) {
            polygon.push_back(polygonLow[low]);
        }
    }

    for (int ml = 0; ml < pMediumLow; ml++) {
        if (! inVector(polygonMediumLowErase, iPMLE, ml)) {
            polygon.push_back(polygonMediumLow[ml]);
        }
    }

    for (int mh = 0; mh < pMediumHigh; mh++) {
        if (! inVector(polygonMediumHighErased, iPMHE, mh)) {
            polygon.push_back(polygonMediumHigh[mh]);
        }
    }

    for (int high = 0; high < pHigh; high++) {
        if (! inVector(polygonHighErase, iPHE, high)) {
            polygon.push_back(polygonHigh[high]);
        }
    }

    return polygon;
}

bool inVector(vector<int>& vector, int size, int elem) {

    for (int i = 0; i < size; i++) {
        if (vector[i] == elem) {
            return true;
        }
    }

    return false;
}




// This only work for the index graph and index file has to be in ascending order
/*
+
X                           XX                 Overload  X
|X Underload               X  X                         X
| X                       X    X  Fitting-             X
|  X                     X      X                     X
|   X                   X        X                   X
|    X                 X          X                 X
|     X    Fitting+   X            X               X
|      X             X     ---------X             X
|-------X           X                X           X
|        X         X                  X         X
|         X       X--------            X       X---------
|          X     X                      X     X
|           X   X                        X   X
|            X X                          X X
--------------X----------------------------X------------+
0           33            50            66           100
*/

vector<Point2D> definePolygonMWL(vector<double>& maxIndexDegree, vector<fuzzy_t>& index) {

    vector<Point2D> polygon;
    int n = maxIndexDegree.size() - 1;

    // Underload
    if (maxIndexDegree[0] > EPSILON) {
        polygon.push_back(Point2D {0, 0});

        if (maxIndexDegree[0] < 1 - EPSILON ) {
            polygon.push_back(Point2D {0, maxIndexDegree[0]});
            polygon.push_back(Point2D {lineX (index[0].p1, index[0].p2, maxIndexDegree[0]),
                                    maxIndexDegree[0]});
            polygon.push_back(index[0].p2);
        } else {
            polygon.push_back(index[0].p1);
            polygon.push_back(index[0].p2);
        }
    }

    // FittingMinus
    if (maxIndexDegree[1] > EPSILON) {

        if (maxIndexDegree[0] < EPSILON) {
            polygon.push_back(index[1].p1);
        }

        if (maxIndexDegree[1] < 1 - EPSILON) {
            polygon.push_back(Point2D{lineX(index[1].p1, index[1].p2, maxIndexDegree[1]), maxIndexDegree[1]});
            polygon.push_back(Point2D{index[1].p2.x, maxIndexDegree[1]});
        } else {
            polygon.push_back(index[1].p2);
        }

        if (maxIndexDegree[2] < EPSILON) {
            polygon.push_back(Point2D{index[1].p2.x, 0});
        }
    }

    // FittingPlus
    if (maxIndexDegree[2] > EPSILON) {

        if (maxIndexDegree[1] < EPSILON) {
            polygon.push_back(Point2D{index[2].p1.x, 0});
        }

        if (maxIndexDegree[2] < 1 - EPSILON) {
            polygon.push_back(Point2D{index[2].p1.x, maxIndexDegree[2]});
            polygon.push_back(Point2D{lineX(index[2].p1, index[2].p2, maxIndexDegree[2]), maxIndexDegree[2]});
        } else {
            if (maxIndexDegree[1] < 1 - EPSILON) {
                polygon.push_back(index[2].p1);
            }
        }

        polygon.push_back(index[2].p2);
    }

    // Overload
    if (maxIndexDegree[3] > EPSILON) {

        if (maxIndexDegree[2] < EPSILON) {
            polygon.push_back(index[3].p1);
        }

        if (maxIndexDegree[3] < 1 - EPSILON) {
            polygon.push_back(Point2D{lineX(index[3].p1, index[3].p2, maxIndexDegree[3]), maxIndexDegree[3]});
            polygon.push_back(Point2D{index[3].p2.x, maxIndexDegree[3]});
        } else {
            polygon.push_back(index[3].p2);
        }

        polygon.push_back(Point2D{index[3].p2.x, 0});
    }

    return polygon;
}

/*
+
|
| Underload                       overload
+-----X             XX           X--------
|      X   medium  X  X medium  X
|       X  lower  X    X upper X
|        X       X      X     X
|         X     X        X   X
|          X   X          X X
|           X X            X
|            X            X X
|           X X          X   X
+----------X---X--------X-----X------------
0       16.5   33     67      83.5
*/

vector<Point2D> definePolygonMWL_HCI(vector<double>& maxIndexDegree,
                                     vector<fuzzy_t>& index) {

    vector<Point2D> polygon;
    int n = maxIndexDegree.size() - 1;

    // UnderLoad: 2Line
    if (maxIndexDegree[0] > EPSILON) {
        polygon.push_back(Point2D {0, 0});

        if (maxIndexDegree[0] < 1 - EPSILON ) {
            polygon.push_back(Point2D {0, maxIndexDegree[0]});
            polygon.push_back(Point2D {lineX (index[0].p2, index[0].p3, maxIndexDegree[0]),
                                    maxIndexDegree[0]});
            polygon.push_back(index[0].p3);
        } else {
            polygon.push_back(index[0].p1);
            polygon.push_back(index[0].p2);
            polygon.push_back(index[0].p3);
        }
    }

    // MediumLower: Line
    if (maxIndexDegree[1] > EPSILON) {

        polygon.push_back(index[1].p1);

        if (maxIndexDegree[1] < 1 - EPSILON) {
            polygon.push_back(Point2D{lineX(index[1].p1, index[1].p2, maxIndexDegree[1]), maxIndexDegree[1]});
            polygon.push_back(Point2D{index[1].p2.x, maxIndexDegree[1]});
        } else {
            polygon.push_back(index[1].p2);
        }

        if (maxIndexDegree[2] < EPSILON) {
            polygon.push_back(Point2D{index[1].p2.x, 0});
        }
    }

    // MediumUpper: Line
    if (maxIndexDegree[2] > EPSILON) {

        if (maxIndexDegree[1] < EPSILON) {
            polygon.push_back(Point2D{index[2].p1.x, 0});
        }

        if (maxIndexDegree[2] < 1 - EPSILON) {
            polygon.push_back(Point2D{index[2].p1.x, maxIndexDegree[2]});
            polygon.push_back(Point2D{lineX(index[2].p1, index[2].p2, maxIndexDegree[2]), maxIndexDegree[2]});
        } else {
            if (maxIndexDegree[1] < 1 - EPSILON) {
                polygon.push_back(index[2].p1);
            }
        }

        polygon.push_back(index[2].p2);
    }

    // Overload
    if (maxIndexDegree[3] > EPSILON) {
        polygon.push_back(index[3].p1);

        if (maxIndexDegree[3] < 1 - EPSILON) {
            polygon.push_back(Point2D{lineX(index[3].p1, index[3].p2, maxIndexDegree[3]), maxIndexDegree[3]});
            polygon.push_back(Point2D{index[3].p3.x, maxIndexDegree[3]});
        } else {
            polygon.push_back(index[3].p2);
            polygon.push_back(index[3].p3);
        }

        polygon.push_back(Point2D{index[3].p3.x, 0});
    }

    // Solve intersections
    for (int i = 0; i < polygon.size() - 2; i++) {
        if (polygon[i].x > polygon[i + 1].x && polygon[i].y < 0.0000001 && polygon[i + 1].y < 0.000001 ) {

            //cout << "1: " << polygon[i - 1].x << ", " << polygon[i - 1].y << endl;
            //cout << "2: " << polygon[i].x << ", " << polygon[i].y << endl;
            //cout << "3: " << polygon[i + 1].x << ", " << polygon[i + 1].y << endl;
            //cout << "4: " << polygon[i + 2].x << ", " << polygon[i + 2].y << endl;

            double x = intersectionX(polygon[i - 1], polygon[i], polygon[i + 1], polygon [i + 2]);
            double y = intersectionY(polygon[i - 1], polygon[i], polygon[i + 1], polygon [i + 2]);

            //cout << "x aqui:" << x << endl;
            //cout << "y aqui:" << y << endl;

            // Hard to explain this. Imagine if the straight line of underload
            // or the straight line of overload go under fitting minus or
            // fitting plus respectively. But I tested =)
            if (polygon[i + 2]. x >= x && polygon[i - 1].x < x) {

                //cout << "(a)" << endl;
                polygon[i].x = x;
                polygon[i].y = y;

                polygon.erase(polygon.begin() + i + 1);

            } else if (polygon[i + 2].x < x) {
                //cout << "(b)" << endl;
                double x = intersectionX(polygon[i - 1], polygon[i], polygon[i + 2], polygon [i + 3]);
                double y = intersectionY(polygon[i - 1], polygon[i], polygon[i + 2], polygon [i + 3]);

                polygon[i].x = x;
                polygon[i].y = y;

                polygon.erase(polygon.begin() + i + 1);
                polygon.erase(polygon.begin() + i + 1);
            } else {
                //cout << "(c)" << endl;
                double x = intersectionX(polygon[i - 2], polygon[i - 1], polygon[i + 1], polygon [i + 2]);
                double y = intersectionY(polygon[i - 2], polygon[i - 1], polygon[i + 1], polygon [i + 2]);

                polygon[i - 1].x = x;
                polygon[i - 1].y = y;
                //cout << "(C)x: " << x << endl;
                //cout << "(C)y: " << y << endl;

                polygon.erase(polygon.begin() + i);
                polygon.erase(polygon.begin() + i);
            }
        }
    }

    return polygon;
}

/*
+
|
| Underload   mediumLowerLoad optimalLoad mediumUpperLoad  overload
+-----X             XX           X           XX           X--------
|      X           X  X         X X         X  X         X
|       X         X    X       X   X       X    X       X
|        X       X      X     X     X     X      X     X
|         X     X        X   X       X   X        X   X
|          X   X          X X         X X          X X
|           X X            X           X            X
|            X            X X         X X          X X
|           X X          X   X       X   X        X   X
+----------X---X--------X---  X-----X-----X------X-----X----------+
      4    5   6    7   8     9  10 11    12 13 14     15 17      20
*/

vector<Point2D> definePolygonMWL_New(vector<double>& maxIndexDegree, vector<fuzzy_t>& index) {

    vector<Point2D> polygon;
    int n = maxIndexDegree.size() - 1;

    // UnderLoad: 2Line
    if (maxIndexDegree[0] > EPSILON) {
        polygon.push_back(Point2D {0, 0});

        if (maxIndexDegree[0] < 1 - EPSILON ) {
            polygon.push_back(Point2D {0, maxIndexDegree[0]});
            polygon.push_back(Point2D {lineX (index[0].p2, index[0].p3, maxIndexDegree[0]),
                                    maxIndexDegree[0]});
            polygon.push_back(index[0].p3);
        } else {
            polygon.push_back(index[0].p1);
            polygon.push_back(index[0].p2);
            polygon.push_back(index[0].p3);
        }
    }

    // MediumLowerLoad, OptimalLoad, MediumUpperLoad: Triangles
    for (int i = 1; i < n; i++) {
        if (maxIndexDegree[i] > EPSILON) {

            polygon.push_back(index[i].p1);

            if (maxIndexDegree[i] < 1 - EPSILON) {
                Point2D middle = {(index[i].p2.x - index[i].p1.x)/2 + index[i].p1.x, 1};
                polygon.push_back(Point2D{lineX(index[i].p1, middle, maxIndexDegree[i]), maxIndexDegree[i]});
                polygon.push_back(Point2D{lineX(middle, index[i].p2, maxIndexDegree[i]), maxIndexDegree[i]});
                polygon.push_back(index[i].p2);
            } else {
                polygon.push_back(Point2D{(index[i].p2.x - index[i].p1.x)/2 + index[i].p1.x, 1});
                polygon.push_back(index[i].p2);
            }
        }
    }

    // High: 2Line
    if (maxIndexDegree[n] > EPSILON) {

        polygon.push_back(index[n].p1);

        if (maxIndexDegree[n] < 1 - EPSILON) {
            polygon.push_back(Point2D{lineX(index[n].p1, index[n].p2, maxIndexDegree[n]), maxIndexDegree[n]});
            polygon.push_back(Point2D{index[n].p3.x, maxIndexDegree[n]});
        } else {
            polygon.push_back(index[n].p2);
            polygon.push_back(index[n].p3);
        }

        polygon.push_back(Point2D{index[n].p3.x, 0});
    }

    /*cout << "polygon antes" << endl << endl;
    printPolygon(polygon);
    cout << "polygon depois" << endl << endl;*/

    // Solve intersections
    for (int i = 0; i < polygon.size() - 2; i++) {
        if (polygon[i].x > polygon[i + 1].x && polygon[i].y < 0.0000001 && polygon[i + 1].y < 0.000001 ) {

            //cout << "aa: " << polygon[i].x << ", " << i << " - " << polygon[i + 1].x << endl;

            double x = intersectionX(polygon[i - 1], polygon[i], polygon[i + 1], polygon [i + 2]);
            double y = intersectionY(polygon[i - 1], polygon[i], polygon[i + 1], polygon [i + 2]);

            polygon[i].x = x;
            polygon[i].y = y;

            polygon.erase(polygon.begin() + i + 1);
        }
    }

    return polygon;
}

// This only work for the index graph and index file has to be in ascending order 
/*
 * No risk     Low risk     Medium Risk    High Risk   Extremely High
 +                                                        risk
 |X               X             X             X             X
 | X             X X           X X           X X           X
 |  X           X   X         X   X         X   X         X
 |   X         X     X       X     X       X     X       X
 |    X       X       X     X       X     X       X     X
 |     X     X         X   X         X   X         X   X
 |      X   X           X X           X X           X X
 |       X X             X             X             X
 |        X             X X           X X           X X
 |       X X           X   X         X   X         X   X
 |      X   X         X     X       X     X       X     X
 +-----X-----X-------X-------X-----X-------X-----X-------X-------+
0      10    15     35       40   60       65   85       90     100
*/
vector<Point2D> definePolygonBiomarkers(vector<double>& maxIndexDegree, vector<fuzzy_t>& index) {

    vector<Point2D> polygon;
    int n = maxIndexDegree.size() - 1;

    // No_risk
    if (maxIndexDegree[0] > EPSILON) {
        if (maxIndexDegree[0] < 1 - EPSILON) {
            polygon.push_back(Point2D {0, 0});
            polygon.push_back(Point2D {0, maxIndexDegree[0]});
            polygon.push_back(Point2D {lineX (index[0].p1, index[0].p2, maxIndexDegree[0]),
                                    maxIndexDegree[0]});
            polygon.push_back(index[0].p2);
        } else {
            polygon.push_back(Point2D {0, 0});
            polygon.push_back(index[0].p1);
            polygon.push_back(index[0].p2);
        }
    }

    for (int i = 1; i < n; i++) {
        if (maxIndexDegree[i] > EPSILON) {
            // Need to add first point
            //if (maxIndexDegree[i - 1] < 0.001) {
            polygon.push_back(index[i].p1);
            //}

            if (maxIndexDegree[i] < 1 - EPSILON) {
                Point2D middle = {(index[i].p2.x - index[i].p1.x)/2 + index[i].p1.x, 1};
                polygon.push_back(Point2D{lineX(index[i].p1, middle, maxIndexDegree[i]), maxIndexDegree[i]});
                polygon.push_back(Point2D{lineX(middle, index[i].p2, maxIndexDegree[i]), maxIndexDegree[i]});
                polygon.push_back(index[i].p2);
            } else {
                polygon.push_back(Point2D{(index[i].p2.x - index[i].p1.x)/2 + index[i].p1.x, 1});
                polygon.push_back(index[i].p2);
            }
        }
    }

    // Extremely_high_risk
    if (maxIndexDegree[n] > EPSILON) {

        //if (maxIndexDegree[n - 1] < 0.001) {
        polygon.push_back(index[n].p1);
        //}

        if (maxIndexDegree[n] < 1 - EPSILON) {
            polygon.push_back(Point2D{lineX (index[n].p1, index[n].p2, maxIndexDegree[n]), maxIndexDegree[n]});
            polygon.push_back(Point2D{100, maxIndexDegree[n]});
        } else {
            polygon.push_back(Point2D{100, 1});
        }

        polygon.push_back(Point2D{100, 0});
    }

    // Solve intersections
    // First three points are the first triangle
    // Last two points are the last line of the last triangle
    for (int i = 3; i < polygon.size() - 2; i++) {
        //cout << "x: " << polygon[i - 1].x << ", y:" << polygon[i - 1].y << endl;
        //cout << "x: " << polygon[i].x << ", y:" << polygon[i].y << endl << endl;
        if (polygon[i - 1].x > polygon[i].x) {
            double y = intersectionY(polygon[i - 2], polygon [i - 1], polygon[i], polygon [i + 1]);
            //cout << "int: " << y << endl;

            double x = (polygon[i - 1].x - polygon[i].x) / 2 + polygon[i].x;

            polygon[i].x = x;
            polygon[i].y = y;

            polygon.erase(polygon.begin() + i - 1);
        }
    }

    return polygon;
}

// This only work for the index graph and index file has to be in ascending order 
/*
 * Low      Medium low     Medium high     High
 +
 |X               X             X             X
 | X             X X           X X           X
 |  X           X   X         X   X         X
 |   X         X     X       X     X       X
 |    X       X       X     X       X     X
 |     X     X         X   X         X   X
 |      X   X           X X           X X
 |       X X             X             X
 |        X             X X           X X
 |       X X           X   X         X   X
 |      X   X         X     X       X     X
 +-----X-----X-------X-------X-----X-------X---
0     0.2   0.3    0.45      0.55 0.7     0.8  1
*/
vector<Point2D> definePolygonTrust_Linear(vector<double>& maxIndexDegree, vector<fuzzy_t>& index) {

    vector<Point2D> polygon;
    int n = maxIndexDegree.size() - 1;

    // Low trust
    if (maxIndexDegree[0] > EPSILON) {
        if (maxIndexDegree[0] < 1 - EPSILON) {
            polygon.push_back(Point2D {0, 0});
            polygon.push_back(Point2D {0, maxIndexDegree[0]});
            polygon.push_back(Point2D {lineX (index[0].p1, index[0].p2, maxIndexDegree[0]),
                                    maxIndexDegree[0]});
            polygon.push_back(index[0].p2);
        } else {
            polygon.push_back(Point2D {0, 0});
            polygon.push_back(index[0].p1);
            polygon.push_back(index[0].p2);
        }
    }

    // Medium low and medium high
    for (int i = 1; i < n; i++) {
        if (maxIndexDegree[i] > EPSILON) {
            // Need to add first point
            //if (maxIndexDegree[i - 1] < 0.001) {
            polygon.push_back(index[i].p1);
            //}

            if (maxIndexDegree[i] < 1 - EPSILON) {
                Point2D middle = {(index[i].p2.x - index[i].p1.x)/2 + index[i].p1.x, 1};
                polygon.push_back(Point2D{lineX(index[i].p1, middle, maxIndexDegree[i]), maxIndexDegree[i]});
                polygon.push_back(Point2D{lineX(middle, index[i].p2, maxIndexDegree[i]), maxIndexDegree[i]});
                polygon.push_back(index[i].p2);
            } else {
                polygon.push_back(Point2D{(index[i].p2.x - index[i].p1.x)/2 + index[i].p1.x, 1});
                polygon.push_back(index[i].p2);
            }
        }
    }

    // High
    if (maxIndexDegree[n] > EPSILON) {

        //if (maxIndexDegree[n - 1] < 0.001) {
        polygon.push_back(index[n].p1);
        //}

        if (maxIndexDegree[n] < 1 - EPSILON) {
            polygon.push_back(Point2D{lineX (index[n].p1, index[n].p2, maxIndexDegree[n]), maxIndexDegree[n]});
            polygon.push_back(Point2D{1, maxIndexDegree[n]});
        } else {
            polygon.push_back(Point2D{1, 1});
        }

        polygon.push_back(Point2D{1, 0});
    }

    // Solve intersections
    // First three points are the first triangle
    // Last two points are the last line of the last triangle
    for (int i = 3; i < polygon.size() - 2; i++) {
        //cout << "x: " << polygon[i - 1].x << ", y:" << polygon[i - 1].y << endl;
        //cout << "x: " << polygon[i].x << ", y:" << polygon[i].y << endl << endl;
        if (polygon[i - 1].x > polygon[i].x) {
            double y = intersectionY(polygon[i - 2], polygon [i - 1], polygon[i], polygon [i + 1]);
            //cout << "int: " << y << endl;

            double x = (polygon[i - 1].x - polygon[i].x) / 2 + polygon[i].x;

            polygon[i].x = x;
            polygon[i].y = y;

            polygon.erase(polygon.begin() + i - 1);
        }
    }

    return polygon;
}

double aggregateRules(vector<rule_t>& rules, string level) {

    double maxOR = 0;
    double sumOR = 0;
    double productOR = 0;

    //      Zadeh: AND = min, OR = max
    //      Godels: AND = min, OR = max
    //      Lukasiewicz: AND = max(a + b - 1, 0)
    //      Product: AND = a.b, OR = a + b - ab

    bool zadehOrGodels = fuzzy_logic.compare("Zadeh") == 0 || fuzzy_logic.compare("Godels") == 0;
    bool lukasiewicz = fuzzy_logic.compare("Lukasiewicz") == 0;
    bool product = fuzzy_logic.compare("Product") == 0;

    for (int rule_i = 0; rule_i < rules.size(); rule_i++) {

        if (rules[rule_i].fuzzy_index.level.compare(level) != 0) {
            continue;
        }

        if (zadehOrGodels && rules[rule_i].mdegree > maxOR) {
            maxOR = rules[rule_i].mdegree;
        } else if (lukasiewicz) {
            sumOR += rules[rule_i].mdegree;
        } else if (product) {
            productOR = productOR + rules[rule_i].mdegree - productOR * rules[rule_i].mdegree;
        }
    }

    if (zadehOrGodels) {
        return maxOR;
    }

    if (lukasiewicz) {
        if (sumOR > 1) {
            sumOR = 1;
        }

        return sumOR;
    }

    if (product) {
        return productOR;
    }
}

void setFuzzyRisk(string level, vector<fuzzy_t>& index, fuzzy_t& fuzzy_index) {

    for (int i = 0; i < index.size(); i++) {
        if (index[i].level.compare(level) == 0) {
            fuzzy_index = index[i];
            break;
        }

        if (i == index.size() - 1) {
            fuzzy_index.level = "not found";
        }
    }
}

void setFuzzyVariable(string name, string level, vector<fuzzy_t>& fuzzy, fuzzy_t& fuzzy_variable) {

    for (int i = 0; i < fuzzy.size(); i++) {
        if (fuzzy[i].name.compare(name) == 0 && fuzzy[i].level.compare(level) == 0) {
            fuzzy_variable = fuzzy[i];
            break;
        }

        if (i == fuzzy.size() - 1) {
            fuzzy_variable.name = "not found";
        }
    }
}

void addFuzzyToRule(string name, string level, vector<fuzzy_t>& fuzzy, rule_t& rule) {

    for (int i = 0; i < fuzzy.size(); i++) {
        if (fuzzy[i].name.compare(name) == 0 && fuzzy[i].level.compare(level) == 0) {
            rule.fuzzy_variables.push_back(fuzzy[i]);
            break;
        }

        if (i == fuzzy.size() - 1) {

            fuzzy_t notFoundFuzzy;
            notFoundFuzzy.name = "not found";

            rule.fuzzy_variables.push_back(notFoundFuzzy);
        }
    }
}

void printFuzzy(vector<fuzzy_t>& fuzzy) {

    for (int i = 0; i < fuzzy.size(); i++) {

        cout << "Name: " << fuzzy[i].name << endl;
        cout << "Level: " << fuzzy[i].level << endl;
        cout << "Type: " << fuzzy[i].shape << endl;

        if (fuzzy[i].shape.compare("triangle") == 0 || fuzzy[i].shape.compare("line") == 0) {
            cout << "Coordinates: " << fuzzy[i].p1.x << "," << fuzzy[i].p1.y << "," 
                                    << fuzzy[i].p2.x << "," << fuzzy[i].p2.y;
        } else if (fuzzy[i].shape.compare("normal") == 0) {
            cout << "Mean: " << fuzzy[i].mean << endl; 
            cout << "Sd: " << fuzzy[i].sd;
        } else if (fuzzy[i].shape.compare("rightNormal") == 0) {
            cout << "Mean: " << fuzzy[i].mean << endl; 
            cout << "Sd: " << fuzzy[i].sd;
        } else if (fuzzy[i].shape.compare("leftNormal") == 0) {
            cout << "Mean: " << fuzzy[i].mean << endl; 
            cout << "Sd: " << fuzzy[i].sd;
        } else if (fuzzy[i].shape.compare("crisp") == 0) {
            cout << "Crisp active: " << fuzzy[i].crispActive;
        } else if (fuzzy[i].shape.compare("2line") == 0) {
           cout << "Coordinates: " << fuzzy[i].p1.x << "," << fuzzy[i].p1.y << "," 
                                   << fuzzy[i].p2.x << "," << fuzzy[i].p2.y << ","
                                   << fuzzy[i].p3.x << "," << fuzzy[i].p3.y;
        }

        cout << endl << endl;
    }
}

void printRulesDegree(vector<rule_t>& rules) {
    for (int i = 0; i < rules.size(); i++) {
        cout << rules[i].label << ": " << rules[i].mdegree << endl;
    }
}

void printRules(vector<rule_t>& rules, vector<attack_t>& attacks, vector<mitigating_t>& mitigatings) {

    for (int i = 0; i < rules.size(); i++) {

        cout << rules[i].label << " - Premisse: ";
        for (int fuzzy_i = 0; fuzzy_i < rules[i].fuzzy_variables.size(); fuzzy_i++) {
            cout << rules[i].fuzzy_variables[fuzzy_i].name << " is " << rules[i].fuzzy_variables[fuzzy_i].level;

            if (fuzzy_i < rules[i].fuzzy_variables.size() - 1) {
                cout << " AND ";
            } else {
                cout << endl;
            }
        }
        cout << "Conclusion: index is " << rules[i].fuzzy_index.level << endl;
    }

    for (int i = 0; i < attacks.size(); i++) {
        if (attacks[i].rebuttal) {
            cout << attacks[i].source << " <=> " << attacks[i].target << endl;
        } else {
            cout << attacks[i].source << " => " << attacks[i].target << endl;
        }
    }

    for (int i = 0; i < mitigatings.size(); i++) {
        for (int j = 0; j < mitigatings[i].source.size(); j++) {

            cout << mitigatings[i].source[j] << ":" << mitigatings[i].sourceLevel[j];

            if (j < mitigatings[i].source.size() - 1) {
                cout << " && ";
            } else {
                cout << " => " << mitigatings[i].target << endl;
            }
        }
    }
}

void printData(vector<vector<string>>& data) {

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            cout << data[i][j];

            if (j != data[i].size() - 1) {
                cout << ",";
            } else {
                cout << endl;
            }
        }
    }
}

void readData(vector<vector<string>>& data, string file) {

    data.resize(number_of_lines(file));

    ifstream infile(file);
    string line;

    string delimiter = ",";
    size_t pos = 0;
    string token;

    int data_i = 0;

    while (getline(infile, line)) {
        while ((pos = line.find(delimiter)) != string::npos) {
            token = line.substr(0, pos);
            data[data_i].push_back(token);
            line.erase(0, pos + delimiter.length());
        }

        data[data_i].push_back(line);
        data_i++;
    }
}

void readRules(vector<fuzzy_t>& fuzzy_variables,
               vector<fuzzy_t>& index,
               vector<rule_t>& rules,
               vector<attack_t>& attacks,
               vector<mitigating_t>& mitigatings,
               string file) {

    ifstream infile(file);
    string line;

    int fuzzy_i = 0;
    int parameter_i = 0;

    string::size_type sz;   // alias of size_t

    while (getline(infile, line)) {

        if (line.size() < 2) {
            continue;
        }

        string attackArrow = " => ";
        string rebuttalArrow = " <=> ";
        string mitigatingDots = ":";

        //cout << "Line: " << line << endl;

        bool rebuttal = line.find(rebuttalArrow) != string::npos;
        bool attack = line.find(attackArrow) != string::npos;
        bool mitigatingAttack = line.find(mitigatingDots) != string::npos;

        if (rebuttal && mitigatingAttack) {
            cout << "Error! Mitigating attack can not be rebuttal at the same time.";
        }

        if (mitigatingAttack && ! attack) {
            cout << "Error! Mitigating attack not well defined.";
        }

        // Attack or rebuttal from one rule to another rule
        if ((attack || rebuttal) && ! mitigatingAttack) {

            size_t pos = line.find(attackArrow);
            if (rebuttal) {
                pos = line.find(rebuttalArrow);
            }

            string source = line.substr(0, pos);
            string target;
            if (rebuttal) {
                target = line.substr(pos + rebuttalArrow.length(), line.length());
            } else {
                target = line.substr(pos + attackArrow.length(), line.length());
            }

            attacks.push_back(attack_t {source, target, rebuttal});

            // Continue to next line in rule file.
            continue;
        }

        if (mitigatingAttack) {
            size_t pos = line.find(attackArrow);

            string target = line.substr(pos + attackArrow.length(), line.length());

            string fullSource = line.substr(0, pos);

            string connector = " && ";
            bool composed = fullSource.find(connector) != string::npos;

            vector<string> source;
            vector<string> sourceLevel;

            string delimiter = ":";
            string token;

            while ((pos = fullSource.find(connector)) != string::npos) {
                token = fullSource.substr(0, pos);

                source.push_back(token.substr(0, token.find(delimiter)));
                sourceLevel.push_back(token.substr(token.find(delimiter) + delimiter.length(), token.length()));

                fullSource.erase(0, pos + connector.length());
            }

            token = fullSource.substr(0, fullSource.length());

            source.push_back(token.substr(0, token.find(delimiter)));
            sourceLevel.push_back(token.substr(token.find(delimiter) + delimiter.length(), token.length()));

            bool error = false;

            for (int i = 0; i < fuzzy_variables.size(); i++) {
                for (int j = 0; j < source.size(); j++) {
                    if (fuzzy_variables[i].name == source[j]) {
                        i = fuzzy_variables.size();
                        break;
                    }

                    if (i == fuzzy_variables.size() - 1) {
                        cout << "Error! Source attack with no associate fuzzy variable!" << endl ; 
                        error = true;
                    }
                }
            }

            for (int j = 0; j < source.size(); j++) {
                if (sourceLevel[j].compare("any") != 0) {
                    for (int i = 0; i < fuzzy_variables.size(); i++) {
                        if (fuzzy_variables[i].name == source[j] && fuzzy_variables[i].level == sourceLevel[j]) {
                            break;
                        }

                        if (i == fuzzy_variables.size() - 1) {
                            cout << "Error! Source level attack with no associate fuzzy variable!" << endl;

                            cout << "Source: " << source[j] << endl;
                            cout << "Source leve: " << sourceLevel[j] << endl;
                            for (int k = 0; k < fuzzy_variables.size(); k++) {
                                cout << "Fuzzy: " << fuzzy_variables[k].name << endl;
                                cout << "Fuzzy: " << fuzzy_variables[k].level << endl;
                            }

                            error = true;
                        }
                    }
                }
            }

            if (! error) {
                mitigatings.push_back(mitigating_t {source, sourceLevel, target,  composed});
            }

            // Continue to next line of the file
            continue;
        }

        rule_t rule;

        // It is not an attack, proceed to read rule
        string delimiter = ",";
        size_t pos = line.find(delimiter);
        if (pos == string::npos) {
            cout << "Error! No conclusion!" << endl ; 
            //exit (0);
        }

        // Save label
        string label = line.substr(0, pos);
        line.erase(0, pos + delimiter.length() + 1);

        // Find comma separating premisses and conclusion
        pos = line.find(delimiter);
        string premisses = line.substr(0, pos);
        // + 1 to skip space
        string conclusion = line.substr(pos + delimiter.length() + 1, line.length());

        // Break premisses
        delimiter = " %is% ";
        size_t n = stringOcurrences(premisses, delimiter);

        if (n == 1) {
            pos = premisses.find(delimiter);
            string variable = premisses.substr(0, pos);
            string level = premisses.substr(pos + delimiter.length(), premisses.length());
            addFuzzyToRule(variable, level, fuzzy_variables, rule);
        } else {

            string connector = " && ";

            vector<fuzzy_t> premissesVariables(stringOcurrences(premisses, connector) + 1);

            size_t pos = 0;
            string token;

            while ((pos = premisses.find(connector)) != string::npos) {
                token = premisses.substr(0, pos);

                string variable = token.substr(0, token.find(delimiter));
                string level = token.substr(token.find(delimiter) + delimiter.length(), token.length());

                addFuzzyToRule(variable, level, fuzzy_variables, rule);

                premisses.erase(0, pos + connector.length());
            }

            token = premisses.substr(0, premisses.length());

            string variable = token.substr(0, token.find(delimiter));
            string level = token.substr(token.find(delimiter) + delimiter.length(), token.length());

            addFuzzyToRule(variable, level, fuzzy_variables, rule);
        }

        pos = conclusion.find(delimiter);

        // *** Only work with indexs
        string level = conclusion.substr(pos + delimiter.length(), conclusion.length());
        setFuzzyRisk(level, index, rule.fuzzy_index);
        rule.label = label;

        bool notFound = false;
        for (int fuzzy_i = 0; fuzzy_i < rule.fuzzy_variables.size(); fuzzy_i++) {
            if (rule.fuzzy_variables[fuzzy_i].name.compare("not found") == 0) {
                cout << "No associated premisse for rule \"" << line << "\"" << endl; 
                notFound = true;
            }
        }

        if (rule.fuzzy_index.name.compare("not found") == 0) {
            cout << "No associated conclusion for rule \"" << line << "\"" << endl;
            notFound = true;
        }

        if (! notFound) {
            rules.push_back(rule);
        }
    }

    for (int attacks_i = 0; attacks_i < attacks.size(); attacks_i++) {

        bool sourceLabel = false;
        bool targetLabel = false;

        for (int rules_i = 0; rules_i < rules.size(); rules_i++) {
            if (attacks[attacks_i].source == rules[rules_i].label) {
                sourceLabel = true;
            }

            if (attacks[attacks_i].target == rules[rules_i].label) {
                targetLabel = true;
            }

            if (! sourceLabel && ! targetLabel && rules_i == rules.size() - 1) {
                cout << "Error! Attack with no associated rule!." << endl;
                cout << "Source: " << attacks[attacks_i].source << endl;
                cout << "Target: " << attacks[attacks_i].target << endl;
            }
        }
    }

    for (int mitigating_i = 0; mitigating_i < mitigatings.size(); mitigating_i++) {

        bool targetLabel = false;

        for (int rules_i = 0; rules_i < rules.size(); rules_i++) {
            if (mitigatings[mitigating_i].target == rules[rules_i].label) {
                targetLabel = true;
            }

            if (! targetLabel && rules_i == rules.size() - 1) {
                cout << "Error! Mitigating attack with no associated rule for target!." << endl;
                cout << "Target: " << mitigatings[mitigating_i].target << endl;
            }
        }
    }
}

void readFuzzy(vector<fuzzy_t>& fuzzy, string file) {

    ifstream infile(file);
    string line;

    int parameter_i = 0;

    /*getline(infile, line);
    string::size_type sz;   // alias of size_t

    if (fuzzy.size() != stoi(line, &sz)) {
        cout << "Error! Too many fuzzy variables!" << endl ;
        cout << "Size of fuzzy: " << fuzzy.size() << endl;
        cout << "Size of file: " << stoi(line, &sz) << endl;
        //exit (0);
    }

    // Skip blank line
    getline(infile, line);*/

    fuzzy_t newFuzzy;
    string::size_type sz;

    while (getline(infile, line)) {

        if (parameter_i == 0) {
            newFuzzy.name = line;
        } else if (parameter_i == 1) {
            newFuzzy.level = line;
        } else if (parameter_i == 2) {
            newFuzzy.shape = line;
        } else if (parameter_i == 3) {

            string delimiter = ",";
            size_t pos = 0;
            string token;
            int coordinate = 0;

            while ((pos = line.find(delimiter)) != string::npos) {
                token = line.substr(0, pos);

                if (newFuzzy.shape.compare("line") == 0 ||
                    newFuzzy.shape.compare("triangle") == 0 ||
                    newFuzzy.shape.compare("2line") == 0) {

                    if (coordinate == 0) {
                        newFuzzy.p1.x = stof(token, &sz);
                    } else if (coordinate == 1) {
                        newFuzzy.p1.y = stof(token, &sz);
                    } else if (coordinate == 2) {
                        newFuzzy.p2.x = stof(token, &sz);
                    } else if (coordinate == 3) {
                        newFuzzy.p2.y = stof(token, &sz);
                    } else if (coordinate == 4) {
                        newFuzzy.p3.x = stof(token, &sz);
                    }
                } else if (newFuzzy.shape.compare("normal") == 0 || 
                           newFuzzy.shape.compare("leftNormal") == 0 ||
                           newFuzzy.shape.compare("rightNormal") == 0) {
                    newFuzzy.mean = stof(token, &sz);
                } // No crisp here because crisps don't have comma

                line.erase(0, pos + delimiter.length());
                coordinate++;
            }

            // Read last coordinate or sd. `While` stops before
            if (newFuzzy.shape.compare("line") == 0 || newFuzzy.shape.compare("triangle") == 0) {
                newFuzzy.p2.y = stof(line, &sz);
            } else if (newFuzzy.shape.compare("normal") == 0 ||
                       newFuzzy.shape.compare("leftNormal") == 0 ||
                       newFuzzy.shape.compare("rightNormal") == 0) {
                newFuzzy.sd = stof(line, &sz);
            } else if (newFuzzy.shape.compare("2line") == 0) {
                newFuzzy.p3.y = stof(line, &sz);
            } else if (newFuzzy.shape.compare("crisp") == 0) {
                newFuzzy.crispActive = stof(line, &sz);
            }

        } else if (parameter_i == 4) {
            parameter_i = 0;
            fuzzy.push_back(newFuzzy);
            newFuzzy = {};
            //fuzzy_i++;
            continue;
        } else {
            cout << "Error! Too many parameters!";
            //exit (0);
        }

        parameter_i++;
    }
}

void setDegree(double input, fuzzy_t& fuzzy_variable) {

    fuzzy_variable.mdegree = 0;

    if (fuzzy_variable.shape.compare("triangle") == 0) {
        double middleX = (fuzzy_variable.p2.x - fuzzy_variable.p1.x) / 2 + fuzzy_variable.p1.x;
        if (input < middleX) {
            fuzzy_variable.mdegree = max(lineY(fuzzy_variable.p1, Point2D {middleX, 1}, input), 0.0);
        } else {
            fuzzy_variable.mdegree = max(lineY(Point2D {middleX, 1}, fuzzy_variable.p2, input), 0.0);
        }
    } else if (fuzzy_variable.shape.compare("2line") == 0) {
        if (input < fuzzy_variable.p1.x) {
            fuzzy_variable.mdegree = 0;
        } else if (input < fuzzy_variable.p2.x) {
            fuzzy_variable.mdegree = lineY(fuzzy_variable.p1, fuzzy_variable.p2, input);
        } else {
            fuzzy_variable.mdegree = min(lineY(fuzzy_variable.p2, fuzzy_variable.p3, input), 1.0);
        }
    } else if (fuzzy_variable.shape.compare("line") == 0) {
        fuzzy_variable.mdegree = lineY(fuzzy_variable.p1, fuzzy_variable.p2, input);
    } else if (fuzzy_variable.shape.compare("normal") == 0) {
        fuzzy_variable.mdegree = max(normalY(fuzzy_variable.mean, fuzzy_variable.sd, input), 0.0);
    } else if (fuzzy_variable.shape.compare("leftNormal") == 0) {
        if (input > fuzzy_variable.mean) {
            fuzzy_variable.mdegree = 0;
        } else {
            fuzzy_variable.mdegree = max(normalY(fuzzy_variable.mean, fuzzy_variable.sd, input), 0.0);
        }
    } else if (fuzzy_variable.shape.compare("rightNormal") == 0) {
        if (input < fuzzy_variable.mean) {
            fuzzy_variable.mdegree = 0;
        } else {
            fuzzy_variable.mdegree = max(normalY(fuzzy_variable.mean, fuzzy_variable.sd, input), 0.0);
        }
    } else if (fuzzy_variable.shape.compare("crisp") == 0) {
        if (input == fuzzy_variable.crispActive) {
            fuzzy_variable.mdegree = 1;
        } else {
            fuzzy_variable.mdegree = 0;
        }
    } else {
        cout << "Error! Variable with wrong type!";
    }

    // Outside function domain
    if (fuzzy_variable.mdegree > 1 || fuzzy_variable.mdegree < 0) {
        fuzzy_variable.mdegree = 0;
    }

    //cout << "Degree = " << fuzzy_variable.mdegree << "\n";
}

void solveRebuttals(vector<rule_t>& rules, attack_t& current) {

    if (! current.rebuttal) {
        return;
    }
    
    int targetIndex = -1;
    int sourceIndex = -1;

    getTargetSourceRuleIndex(targetIndex, sourceIndex, current, rules);

    double sourceDegree = 1 - rules[sourceIndex].rebuttalMDegree;
    double targetDegree = rules[targetIndex].mdegree;

    if (sourceDegree < targetDegree) {
        rules[targetIndex].mdegree = sourceDegree;
    }

    // Apply attack with swapped indexes
    sourceDegree = 1 - rules[targetIndex].rebuttalMDegree;
    targetDegree = rules[sourceIndex].mdegree;

    if (sourceDegree < targetDegree) {
        rules[sourceIndex].mdegree = sourceDegree;
    }
}

void findLeaves(attack_t& current,
                int currentIndex,
                vector<attack_t>& attacks,
                vector<rule_t>& rules,
                vector<int>& leaves) {

    if (current.rebuttal) {
        return;
    }

    for (int i = 0; i < attacks.size(); i++) {
        // Find an attack whose source is the target of the current attack
        // and it is not the same attack.
        if (current.target == attacks[i].source && i != currentIndex && ! attacks[i].rebuttal) {
            findLeaves(attacks[i], i, attacks, rules, leaves);
            break;
        }

        if (i == attacks.size() - 1) {
            //cout << "Index: " << currentIndex << endl;
            //cout << "Source: " << current.source << " - " << current.target;
            leaves.push_back(currentIndex);
        }
    }
}

void getTargetSourceRuleIndex(int& targetIndex,
                              int& sourceIndex,
                              attack_t& attack,
                              vector<rule_t>& rules) {

    // Find source and target rules
    targetIndex = -1;
    sourceIndex = -1;

    for (int i = 0; i < rules.size(); i++) {
        if (rules[i].label == attack.target) {
            targetIndex = i;
        }

        if (rules[i].label == attack.source) {
            sourceIndex = i;
        }
    }

    if (sourceIndex == -1 || targetIndex == -1) {
        cout << "Error! Attack target or source with no associated rule!" << endl;
        cout << "Source: (" << attack.source << ")" << endl;
        cout << "Target: (" << attack.target << ")" << endl;
    }
}

void solveAttacks(attack_t& current,
                  vector<rule_t>& rules,
                  vector<attack_t>& attacks,
                  vector<fuzzy_t>& fuzzy_variables) {

    if (current.rebuttal) {
        return;
    }

    // Check if source is target. If it is need to do recursive call
    for (int i = 0; i < attacks.size(); i++) {
        if (attacks[i].target == current.source && attacks[i].source != current.target && ! attacks[i].rebuttal) {
            solveAttacks(attacks[i], rules, attacks, fuzzy_variables);
        }
    }

    int targetIndex = -1;
    int sourceIndex = -1;

    getTargetSourceRuleIndex(targetIndex, sourceIndex, current, rules);

    double attackerDegree = 1 - rules[sourceIndex].mdegree;

    if (attackerDegree < rules[targetIndex].mdegree) {
        rules[targetIndex].mdegree = attackerDegree;
    }
}

void solveMitigatings(vector<rule_t>& rules,
                      mitigating_t& mitigating,
                      vector<fuzzy_t>& fuzzy_variables) {

    int targetIndex;

    for (int i = 0; i < rules.size(); i++) {
        if (rules[i].label == mitigating.target) {
            targetIndex = i;
            break;
        }

        if (i == rules.size() - 1) {
            cout << "Error! Mitigating target with no associated rule!" << endl;
            cout << "Target: " << mitigating.target;
        }
    }

    if (rules[targetIndex].mdegree < EPSILON) {
        return;
    }

    double degreeSource;

    bool zadehOrGodels = fuzzy_logic.compare("Zadeh") == 0 || fuzzy_logic.compare("Godels") == 0;
    bool lukasiewicz = fuzzy_logic.compare("Lukasiewicz") == 0;
    bool product = fuzzy_logic.compare("Product") == 0;

    if (zadehOrGodels) {
        degreeSource = 1;
    }

    if (lukasiewicz) {
        degreeSource = 0;
    }

    for (int j = 0; j < mitigating.source.size(); j++) {
        for (int k = 0; k < fuzzy_variables.size(); k++) {
            if (fuzzy_variables[k].name.compare(mitigating.source[j]) != 0 ||
                fuzzy_variables[k].level.compare(mitigating.sourceLevel[j]) != 0) {
                continue;
            }

            if (zadehOrGodels) {
                // Taking the min cause it is assumed that 
                // premissed are connected only by AND
                if (fuzzy_variables[k].mdegree < degreeSource) {
                    degreeSource = fuzzy_variables[k].mdegree;
                }
            } else if (lukasiewicz) {
                degreeSource += fuzzy_variables[k].mdegree;
            } else if (product) {
                if (j == 0) {
                    degreeSource = fuzzy_variables[k].mdegree;
                } else {
                    degreeSource *= fuzzy_variables[k].mdegree;
                }
            }
        }
    }

    if (lukasiewicz) {
        degreeSource = degreeSource - (mitigating.source.size() - 1);

        if (degreeSource < 0) {
            degreeSource = 0;
        }
    }

    double attackerDegree = 1 - degreeSource;

    if (attackerDegree < rules[targetIndex].mdegree) {
        rules[targetIndex].mdegree = attackerDegree;
    }
}

double getValue(vector<string>& data, vector<string>& header, string name) {

    string::size_type sz;
    //cout << "\n" << name << "\n";
    for (int i = 0; i < header.size(); i++) {
        //cout << "\n" << header << "\n";
        if (header[i].compare(name) == 0) {

            if (data[i].compare("?") == 0) {
                return -999;
            }

            return stof(data[i], &sz);
        }
    }

    return -999;
}

double getWeight(vector<string>& data, vector<string>& header, string name) {

    string::size_type sz;
    for (int i = 0; i < header.size(); i++) {

        bool isWeight = header[i].find("Weight_") != string::npos && 
                        header[i].find(name) != string::npos;

        if (isWeight) {

            if (data[i].compare("?") == 0) {
                return 1;
            }
            return stof(data[i], &sz);
        }
    }

    return 1;
}

void normalizeWeights(vector<fuzzy_t>& fuzzy_variables) {

    double min = fuzzy_variables[0].weight;
    double max = fuzzy_variables[0].weight;

    for (int i = 0; i < fuzzy_variables.size(); i++) {

        if (fuzzy_variables[i].weight > max) {
            max = fuzzy_variables[i].weight;
        }

        if (fuzzy_variables[i].weight < min) {
            min = fuzzy_variables[i].weight;
        }
    }

    if (max - min == 0) {
        for (int i = 0; i < fuzzy_variables.size(); i++) {
            fuzzy_variables[i].weight = 1;
        }

        return;
    }

    for (int i = 0; i < fuzzy_variables.size(); i++) {
        fuzzy_variables[i].weight = (fuzzy_variables[i].weight - min) / (max - min);
    }
}

void applyWeights(vector<rule_t>& rules) {

    for (int rule_i = 0; rule_i < rules.size(); rule_i++) {

        double maxWeight = rules[rule_i].fuzzy_variables[0].weight;
        for (int fuzzy_i = 1; fuzzy_i < rules[rule_i].fuzzy_variables.size(); fuzzy_i++) {
            if (rules[rule_i].fuzzy_variables[fuzzy_i].weight > maxWeight) {
                maxWeight = rules[rule_i].fuzzy_variables[fuzzy_i].weight;
            }
        }

        rules[rule_i].mdegree *= maxWeight;
    }
}


// Given two points, define a line and finds y correspondent to x
double lineY (Point2D p1, Point2D p2, double x) {

    double a = (p2.y - p1.y)/(p2.x - p1.x);

    double b = -(p1.x * a) + p1.y;

    return (a*x + b);
}

// Given two points, define a line and finds x correspondent to y
double lineX (Point2D p1, Point2D p2, double y) {

    double a = (p2.y - p1.y)/(p2.x - p1.x);

    double b = -(p1.x * a) + p1.y;

    return (y - b)/a;
}


// https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
double intersectionY(Point2D p1, Point2D p2, Point2D p3, Point2D p4) {

    double numerator = (p1.x * p2.y - p1.y * p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x * p4.y - p3.y * p4.x);

    double denominador = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);

    return numerator / denominador;
}

// https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
double intersectionX(Point2D p1, Point2D p2, Point2D p3, Point2D p4) {

    double numerator = (p1.x * p2.y - p1.y * p2.x) * (p3.x - p4.x) - (p1.x - p2.x) * (p3.x * p4.y - p3.y * p4.x);

    double denominador = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);

    return numerator / denominador;
}


double normalY(double mean, double sd, double x) {
    // Last part "((1 / (sqrt(2 * 3.14159) / sd)))" is to normalize height to [0, 1]
    // First part normal distribution formula http://seankross.com/notes/dpqr/
    return ((1 / (sqrt(2 * 3.14159) / sd)) * pow(2.71828, - pow(x - mean, 2)/(2 * sd * sd))) / ((1 / (sqrt(2 * 3.14159) / sd)));
}

double normalXLeft(double mean, double sd, double y) {

    //cout << "Mean: " << mean << endl;
    //cout << "Sd: " << sd << endl;
    //cout << "y: " << y << endl;
   double test;
   for (double x = mean -4 * sd; x < mean; x += 0.01) {
        //cout << "Normal y: " << normalY(mean, sd, x) << endl;
        if (abs(y - normalY(mean, sd, x)) < 0.01) {
            test = x;
        }
    }

    double delta = sqrt(4 * mean * mean - 4 * (mean * mean + log(y) * 2 * sd * sd));

    double root1 = (2 * mean + delta) / 2;
    double root2 = (2 * mean - delta) / 2;

    //cout << "r1: " << root1 << endl;
    //cout << "r2: " << root2 << endl;
    //cout << "test: " << test << endl;

    if (abs(y - normalY(mean, sd, root1)) < 0.01 && root1 <= mean) {
        return root1;
    }

    if (abs(y - normalY(mean, sd, root2)) < 0.01 && root2 <= mean) {
        return root2;
    }

    cout << "Error! Couldn't find normalXLeft!" << endl;
    return -999;
}

double normalXRight(double mean, double sd, double y, double precision) {

    //cout << "Mean: " << mean<< endl;
    //cout << "Sd: " << sd<< endl;
    //cout << "Top: " << mean +4 * sd<< endl;
    //cout << "y: " << y<< endl;

    /*double x = mean;
    do {

        if (abs(y - normalY(mean, sd, x)) < 0.01) {
            return x;
        }

        if (y - normalY(mean, sd, x) > 0.01) {
            x += 100 * precision;

            if (x >= mean +4 * sd) {
                x = mean +4 * sd;
            }
        }

        if (y - normalY(mean, sd, x) <= 0.01) {
            x -= precision;
        }

    } while (abs(y - normalY(mean, sd, x)) > 0.01);*/

    /*double test;
    for (double x = mean; x < mean +4 * sd; x += precision) {
        //cout << "Normal y: " << normalY(mean, sd, x) << endl;
        if (abs(y - normalY(mean, sd, x)) < 0.01) {
            //return x;
            test = x;
        }
    }*/

    double delta = sqrt(4 * mean * mean - 4 * (mean * mean + log(y) * 2 * sd * sd));

    double root1 = (2 * mean + delta) / 2;
    double root2 = (2 * mean - delta) / 2;

    /*cout << "r1: " << root1 << endl;
    cout << "r2: " << root2 << endl;
    cout << "test: " << test << endl;*/

    if (abs(y - normalY(mean, sd, root1)) < 0.01 && root1 >= mean) {
        return root1;
    }

    if (abs(y - normalY(mean, sd, root2)) < 0.01 && root2 >= mean) {
        return root2;
    }

    cout << "Error! Couldn't find normalXRight!" << endl;
    return -999;
}


double phi (double x) {
    return ((1/sqrt(2*M_PI))*exp(-0.5*x*x));
}

// https://stackoverflow.com/questions/2792443/finding-the-centroid-of-a-polygon
Point2D compute2DPolygonCentroid(vector<Point2D>& vertices) {
    Point2D centroid = {0, 0};
    double signedArea = 0.0;
    double x0 = 0.0; // Current vertex X
    double y0 = 0.0; // Current vertex Y
    double x1 = 0.0; // Next vertex X
    double y1 = 0.0; // Next vertex Y
    double a = 0.0;  // Partial signed area

    // For all vertices except last
    int i=0;
    for (i=0; i<vertices.size()-1; ++i)
    {
        x0 = vertices[i].x;
        y0 = vertices[i].y;
        x1 = vertices[i+1].x;
        y1 = vertices[i+1].y;
        a = x0*y1 - x1*y0;
        signedArea += a;
        centroid.x += (x0 + x1)*a;
        centroid.y += (y0 + y1)*a;
    }

    // Do last vertex separately to avoid performing an expensive
    // modulus operation in each iteration.
    x0 = vertices[i].x;
    y0 = vertices[i].y;
    x1 = vertices[0].x;
    y1 = vertices[0].y;
    a = x0*y1 - x1*y0;
    signedArea += a;
    centroid.x += (x0 + x1)*a;
    centroid.y += (y0 + y1)*a;

    signedArea *= 0.5;
    centroid.x /= (6.0*signedArea);
    centroid.y /= (6.0*signedArea);

    return centroid;
}

int stringOcurrences(string s, string target) {
   int occurrences = 0;
   string::size_type pos = 0;
   while ((pos = s.find(target, pos )) != string::npos) {
          ++ occurrences;
          pos += target.length();
   }

   return occurrences;
}


int number_of_lines(string file) {

    int result = 0;

    std::string line;
    std::ifstream myfile(file);

    while (std::getline(myfile, line))
        ++result;

    return result;
}
