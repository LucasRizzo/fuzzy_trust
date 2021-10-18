## Non-monotonic fuzzy reasoning implementation

Custom implementation of fuzzy reasoning models, employing a Mandani fuzzy inference process [[1]](#1), and a non-monotonic layer through Possibility Theory [[2]](#2) for resolving contractions. While there are other better developed and documented fuzzy packages, such as FuzzyR in the R language, no other implementation with a non-monotonic layer was previously found till the creation of this repository. The example provided is related to an application of computational trust, but this implementation has been used in many other publications also focused on different applications [[3]](#3),[[4]](#4),[[5]](#5),[[6]](#6).

#### Fuzzy membership functions of rule's variables

The following shapes have been defined: triangle, two lines, line, normal, right normal (normal without left part), and left normal (normal without the right part). Normal shapes are implemented with 1000 points. More points can be added but this can add considerably the execution time during the defuzzification stage.
Text documents are available to exemplify the use of such shapes.

#### T-norms and T-Conorms

The Zadeh, Lukasiewicz and Product fuzzy operators have been implemented, including the associativity for more than two variables. This can be set in the global variables in the .cpp file.

#### Fuzzy membership functions of rule's consequents

Two membership functions have been hardcoded for rule's consequents based on the illustrations below. Other membership functions for rule's consequents also need to be hardcoded in order to find the outer envelope of the truncated membership functions.

![Screenshot from 2021-07-09 13-08-41](https://user-images.githubusercontent.com/14744131/125075678-d9579f00-e0b6-11eb-9ce9-3f3273cfa2db.png)
![Screenshot from 2021-07-09 13-08-33](https://user-images.githubusercontent.com/14744131/125075683-dbb9f900-e0b6-11eb-960d-6eddd08ff240.png)

#### Defuzzification methods

Two methods have been implemented for defuzzification: centroid and mean of max. Other methods also need to be hardcoded.


## How to use

The following files are required:

##### A parameters file
An example can be seen in the ```parameters.txt``` file. This will be used to define the membership function of fuzzy variables levels that can be employed by fuzzy rules.
##### A consequent file
An example can be seen in the ```trust_index.txt``` file. This follows the same syntax of the parameters file, but is used to define only the fuzzy membership functions of the fuzzy variables applied as the consequent of rules.
##### A rules file
Examples can be seen in the ```rules.txt``` file. This will be used to define the fuzzy rules and contradictions. An example of rule is:
```
R1, age %is% age.high, risk %is% risk.high
```
In R1, a fuzzy variable ```age``` infers ```risk``` at level ```risk.high``` if ```age``` is level ```age.high```. ```age.high``` and ```risk.high``` need to be defined in the parameters file. A logical ```AND``` operator is implemented to add more than one fuzzy variable in the premises with the double and syntax: ```&&```.
Two rules can contradict each other with the following notation:
```
R1 => R2
```
This implies that if R1 is true then R2 is being contradicted and needs to be reevaluated. A mutual contradiction can also be defined as:
```
R1 <=> R2
```
Lastly, a set of premises can also contradict a rule with the following notation:
```
health:health.low => R1
```
```health.low``` needs to be defined in the parameters file. The logical AND operator ```&&``` can also be used here to add more than one premise.

##### A dataset file
A csv file whose headers match the name of the variables. Each row will be an instance to be tested by the fuzzy reasoning model. Weights for features can be defined using the syntax ```Weight_<fuzzy_variable>``` in the file's header. The file ```portuguese.csv``` gives an example. Headers that don't match any fuzzy variable will be ignored. 

##### A models file
A models file is defined in order to facilitate the execution of models with different configurations. An example can be seen in the file ```models.txt```.
Each line in the file follows the syntax:
```Model_Names,Fuzzy_logic,Weights,Type```
This configuration will return two results: one with the centroid and one with the mean max defuzzification approach (this could be easily changed in the code)
1. Model_names: name of the models and file in which results will be exported.
2. Fuzzy_logic: Zadeh, Lukasiewicz or Product
3. Weights: Yes or No for rule weights.
4. Type: This parameter allows the separation of models in different folders. In the example there was one folder for models using linear membership functions, and another for models using gaussian membership functions. It allows the execution of models with the same configuration, but different files for rules, parameters, and consequents.

## Contact

Contact me at lucasmrizzo@gmail.com


## References
<a id="1">[1]</a> 
E. H. Mamdani, Application of fuzzy algorithms for control of simple
dynamic plant, Proceedings of the Institution of Electrical Engineers 121
(1974) 1585–1588.

<a id="2">[2]</a>
W. Siler, J. J. Buckley, Fuzzy expert systems and fuzzy reasoning, John
Wiley & Sons, 2005

<a id="3">[3]</a>
L. Rizzo, Evaluating the Impact of Defeasible Argumentation as a Modelling Technique for Reasoning under Uncertainty, Ph.D. thesis, Technological University Dublin, 2020.

<a id="4">[4]</a>
L. Rizzo, L. Majnaric, L. Longo,
A comparative study of defeasible
argumentation and non-monotonic fuzzy reasoning for elderly survival
prediction using biomarkers, in: C. Ghidini, B. Magnini, A. Passerini,
P. Traverso (Eds.), AI*IA 2018 – Advances in Artificial Intelligence,
Springer International Publishing, Cham, 2018, pp. 197–209.

<a id="5">[5]</a>
L. Rizzo, L. Longo, An empirical evaluation of the inferential capacity
of defeasible argumentation, non-monotonic fuzzy reasoning and expert
systems, Expert Systems with Applications 147 (2020) (in press).

<a id="6">[6]</a>
Longo, L. and Rizzo, L., Examining the modelling capabilities of defeasible argumentation and non-monotonic fuzzy reasoning Knowledge-Based Systems, (2021) (in press).
