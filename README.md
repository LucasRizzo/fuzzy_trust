## Non-monotonic fuzzy reasoning implementation

Custom implementation of fuzzy reasoning models, employing a Mandani fuzzy inference process [[1]](#1), and a non-monotonic layer through Possibility Theory [[2]](#2) for resolving contractions. While there are other better developed and documented fuzzy packages, such as FuzzyR in the R language, no other implementation with a non-monotonic layer was previously found to the creation of this repository. The example provided related to an application of computational trust, but this implementation has been used in many other publications also focused on different applications [[3]](#3),[[4]](#4),[[5]](#5),[[6]](#6).

#### Fuzzy membership functions of rule's variables

The following shapes have been defined: triangle, two lines, line, normal, right normal (normal without left part), and left normal (normal without the right part).
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

Three text files are required:

##### A parameters file
An example can be seen in the ```parameters.txt``` file. This will be used to define the membership function of fuzzy variables levels that can be employed by fuzzy rules.
##### A consequent file
An example can be seen in the ```trust_index.txt``` file. This follows the same syntax of the parameters file, but is used to define only the fuzzy membership functions of the fuzzy variables applied as the consequent of rules.
##### A rules file
An example can be seen in the ```rules.txt``` file. This will be used to define the fuzzy rules and contradictions. An example is:
```
R1, age %is% age.high, risk %is% risk.high
```
This examples defined a rule R1, in which the fuzzy variable ```age``` infers ```risk``` at level ```risk.high``` if ```age``` is level ```age.high```. ```age.high``` and ```risk.high``` need to be defined in the parameters file. A logical AND operator is implemented to add more than one fuzzy variable in the premisses with the double and syntax: ```&&```.
Two rules can contradict each other with the following notation:
```
R1 => R2
```
This implied that if R1 is true than R2 is being contracted and needs to be reavaluated. A mutual contradiction can also be defined as:
```
R1 <=> R2
```
Lastly, a set of premisses can also contradict a rule with the following notation:
```
health:health.low => R1
```
```health.low``` needs to be defined in the parameters file. The logical AND operator ```&&``` can also be used here to add more than one premisse.

### References
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