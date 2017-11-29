# MetQuest

MetQuest is a dynamic programming based algorithm for identifying all possible pathways from metabolic networks between the source and the target metabolites.

![Image of MetQuest](https://github.com/aarthi31/MetQuest/blob/master/Images/GitHUbMetQuest.png)

## Getting started

1. Clone this repository to your computer using ```git``` or [download the repository](https://github.com/aarthi31/MetQuest/) and decompress it. All the executable codes can be found in the folder Codes.  
2. Install [Python 2.7.14](https://www.python.org/downloads/)

Ensure that the following packages (with the specified version) are installed:
[NetworkX 1.9.1](http://networkx.github.io/), [cobrapy 0.5.4](https://github.com/opencobra/cobrapy) and
[libSBML 5.13.0](http://sbml.org/Software/libSBML/docs/python-api/libsbml-downloading.html)

Use the following commands on Terminal (Mac and Linux) to install the packages
```pip install cobra==0.5.4```

```pip install scipy==0.13.3```

```pip install networkx==1.9.1```

```pip install python-libsbml==5.13.0```

On windows, navigate to folder ```C:\Python27\Scripts\```, type ```cmd``` in ```Win + R```. Type

```pip.exe install cobra==0.5.4```

```pip.exe install scipy==0.13.3```

```pip.exe install networkx==1.9.1```

```pip.exe install python-libsbml==5.13.0```


## Input

Folder whose structure is as shown:
```
   
    ├── Example                 # Folder  
    │   ├── SBML model(s) of metabolic networks          # XML files of the metabolic networks (COBRA-compatible)
    │   ├── seed_mets.txt         # Text file containing the seed metabolites separated by a newline
    │   ├── source_mets.txt       # Text file containing the source metabolites separated by a newline
    │   └── target_mets.txt       # Text file containing the target metabolites separated by a newline
    └── ...
 ```

## Running MetQuest

On Mac and Linux

Go to terminal and navigate to the directory where the codes are downloaded. For instance, if the codes are downloaded in Desktop, go to terminal and type 
``` cd Desktop/MetQuest/Codes/``` 

To run MetQuest, type 

``` python exec_metquest.py ```

Follow the instructions on the screen.

On Windows

Open Python IDE, open the ```exec_metquest.py``` and execute it. Follow the instructions on the screen.

## Output

The output pathways get saved as text files in a new folder called ```Results```. The output in the files follow the format:

Reaction #

Reaction Name: reaction



