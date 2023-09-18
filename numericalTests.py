import json 
import computeSolvent

def make_input(charges, charge_coords, perm_formulation, solvent_type, test_case):

    in_dict = {
        "order" : 9 ,
        "box" : [-16,16],
        "prec" : 1.0e-5,
        "charge_width" : 1000,
        "charges" : charges,
        "charge_coords" : charge_coords,
        "cav_coords" : [[0.0, 0.0, 0.0]],
        "cav_radii" : [3.7794522509156563],  # 2.0 bohr
        "boundary_width" : 0.2,
        "eps_out" : 78.54,
        "perm_formulation" : perm_formulation,
        "solvent_type" : solvent_type,
        "ionic_strength" : 0.1,
        "kain_hist" : 0,
        "test_case" : test_case
    }
    return in_dict


solvers = ["gpe", "pb", "lpb"]
formulations = ["exponential", "linear"]
cases = [
    {
        "charges" : [1.0],
        "charge_coords" : [[0.0, 0.0, 0.0]]
    },
    {
        "charges" : [1.0, 1.0],
        "charge_coords" : [[1.8897261246257702, 0.0, 0.0],
                           [-1.8897261246257702, 0.0, 0.0]]
    },
    {
        "charges" : [1.0, 1.0, -1.0, -1.0],
        "charge_coords" : [[1.8897261246257702, 0.0, 0.0],
                           [-1.8897261246257702, 0.0, 0.0],
                           [0.0, 1.8897261246257702, 0.0],
                           [0.0, -1.8897261246257702, 0.0]]
    },
    {
        "charges" : [1.0, 1.0, -1.0, -1.0],
        "charge_coords" : [[2.267671349550924, 0.0, 0.0],
                           [-2.267671349550924, 0.0, 0.0],
                           [0.0, 2.267671349550924, 0.0],
                           [0.0, -2.267671349550924, 0.0]]
    },
    {
        "charges" : [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        "charge_coords" : [[0.7558904498503081, 0.0, 0.0],
                           [0.0, 1.5117808997006161, 0.0],
                           [0.0, 0.0, 2.267671349550924],
                           [0.0, 0.0, -0.7558904498503081],
                           [-1.5117808997006161, 0.0, 0.0],
                           [0.0, -2.267671349550924, 0.0]]
    },
    {
        "charges" : [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        "charge_coords" : [[0.37794522492515403, 0.37794522492515403, 0.37794522492515403],
                           [0.9448630623128851, 0.9448630623128851, 0.9448630623128851],
                           [1.5117808997006161, 1.5117808997006161, 1.5117808997006161],
                           [-0.37794522492515403, 0.37794522492515403, -0.37794522492515403],
                           [0.9448630623128851, -0.9448630623128851, 0.9448630623128851],
                           [-1.5117808997006161, -1.5117808997006161, -1.5117808997006161]]
    }
]
numeric_tests = {}
print("Running numeric tests")
for solver in solvers:
    print("-Running tests for solver: ", solver)
    numeric_tests[solver]  = {}
    for fomulation in formulations:
        print("--Running tests for formulation: ", fomulation)
        numeric_tests[solver][fomulation] = []
        for i, case in enumerate( cases):
            print("---Running tests for case: ", i)
            print("----creating input")
            in_dict = make_input(**cases[i], perm_formulation=fomulation, solvent_type=solver, test_case=i)
            print("----input created")
            print("----running solver")
            E_r, iterations = computeSolvent.run(**in_dict)
            output = {"E_r": E_r, "iterations": iterations}
            print("----solver finished, output: ", output)
            print("---------------------------------\n")
            numeric_tests[solver][fomulation].append({"input": in_dict, "output": output})

            

with open("numeric_tests.json", "w") as f:
    json.dump(numeric_tests, f, indent=4)
