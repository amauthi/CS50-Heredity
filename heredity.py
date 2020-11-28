import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }
       

    # Loop over all sets of people who might have the trait
    names = set(people)
    
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue


        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):
                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)
                
    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.*
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]

# doesn't work with Harry !! because of its parents ?
def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and*
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    
    print("one_gene = ",one_gene)
    print("two_genes =",two_genes)
    print("have_trait = ",have_trait)
    
    # contains all probas for the different people
    probEveryone = []
    
    for p in people:
        # print("people = ", p)
        # print("person =", people[p])
    
        # for a single person, contains all probas
        probSinglePerson = []
        
        person_gene = 0
        has_person_trait = False
        
        # for all people, let's check in which category they are
        probPerso0Gene = 0
        probPerso1Gene = 0
        probPerso2Gene = 0
        
        
        # if the person has a mother and a father
        if people[p]["mother"] != None and people[p]["father"] != None:
            # for each case :
            # P0 = P(0fromM)*P(0fromF)
            # P1 = P(0fromM)*P(1fromF) + P(1fromM)*P(0fromF)
            # P2 = P(1fromM)*P(1fromF)
            
            # the mother has one gene
            if people[p]["mother"] in one_gene:
                # the father has one gene
                if people[p]["father"] in one_gene:
                    probPerso0Gene = 0.51*0.51
                    probPerso1Gene = 0.51*0.49 + 0.49*0.51
                    probPerso2Gene = 0.49*0.49
                # the father has two genes
                elif people[p]["father"] in two_genes:
                    # 0.01*0.51 because we need 1 mutation and to choose the mutated gene
                    probPerso0Gene = 0.51*PROBS["mutation"]
                    probPerso1Gene = 0.51*0.99 + 0.49*PROBS["mutation"]
                    probPerso2Gene = 0.49*0.99
                # father has no gene
                else:
                    probPerso0Gene = 0.51*0.99
                    probPerso1Gene = 0.51*PROBS["mutation"] + 0.49*0.99
                    probPerso2Gene = 0.49*PROBS["mutation"]
            # the mother has two genes
            elif people[p]["mother"] in two_genes:
                # the father has one gene
                if people[p]["father"] in one_gene:
                    probPerso0Gene = PROBS["mutation"]*0.51
                    probPerso1Gene = PROBS["mutation"]*0.49 + 0.99*0.51
                    probPerso2Gene = 0.99*0.49
                # the father has two genes
                elif people[p]["father"] in two_genes:
                    probPerso0Gene = PROBS["mutation"]*0.51*PROBS["mutation"]*0.51
                    probPerso1Gene = PROBS["mutation"]*0.99 + 0.99*PROBS["mutation"]
                    probPerso2Gene = 0.99*0.99
                # father has no gene   
                else:
                    probPerso0Gene = PROBS["mutation"]*0.99
                    probPerso1Gene = PROBS["mutation"]*0.51 + 0.99*0.99
                    probPerso2Gene = 0.99*PROBS["mutation"]
            # the mother has no gene
            else:
                # the father has one gene
                if people[p]["father"] in one_gene:              
                    probPerso0Gene = 0.99*0.51
                    probPerso1Gene = 0.99*0.49 + PROBS["mutation"]*0.51
                    probPerso2Gene = PROBS["mutation"]*0.49
                # the father has two genes
                elif people[p]["father"] in two_genes:
                    probPerso0Gene = 0.99*PROBS["mutation"]
                    probPerso1Gene = 0.99*0.99 + PROBS["mutation"]*PROBS["mutation"]
                    probPerso2Gene = PROBS["mutation"]*0.99
                # father has no gene
                else:
                    # nobody has genes
                    probPerso0Gene = 0.99*0.99
                    probPerso1Gene = 0.99*PROBS["mutation"] + PROBS["mutation"]*0.99
                    probPerso2Gene = PROBS["mutation"]*PROBS["mutation"]
                    
        
            print("proba of 0 gene for ", p,"  =",probPerso0Gene)
            print("proba of 1 gene for ", p,"  =",probPerso1Gene)
            print("proba of 2 genes for ", p,"  =",probPerso2Gene)
            

        # the person has no parent
        else:
            probPerso0Gene = PROBS["gene"][0]
            probPerso1Gene = PROBS["gene"][1]
            probPerso2Gene = PROBS["gene"][2]
                    
                    
        if p in one_gene:
            person_gene = 1
            probSinglePerson.append(probPerso1Gene)
        elif p in two_genes:
            person_gene = 2
            probSinglePerson.append(probPerso2Gene)
        else:
            person_gene = 0
            probSinglePerson.append(probPerso0Gene)
        
        
        # calculate the prob to have trait or no based on the number of gene
        if p in have_trait:         
            has_person_trait = True
            probSinglePerson.append(PROBS["trait"][person_gene][True])
                
            
        if has_person_trait == False:
            probSinglePerson.append(PROBS["trait"][person_gene][False])
        
        # print("all probas for ", p, " = ",probSinglePerson)
        
        
        resultSinglePerson = 1
        for prob in probSinglePerson:
            resultSinglePerson = prob*resultSinglePerson
            
        
            
        probEveryone.append(resultSinglePerson)
        
    
    # let's multiply all probByPerson
    resultFinal = 1
    for probFinal in probEveryone:
       resultFinal = probFinal*resultFinal
    
    
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! joint proba =", resultFinal)
    return resultFinal


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    
#    print("probabilities before update =",probabilities)
    
    for key in probabilities:
        if key in one_gene:
            probabilities[key]['gene'][1] = probabilities[key]['gene'][1] + p
        elif key in two_genes:
            probabilities[key]['gene'][2] = probabilities[key]['gene'][2] + p
        else:
           probabilities[key]['gene'][0] = probabilities[key]['gene'][0] + p
        
        if key in have_trait:
            probabilities[key]['trait'][True] = probabilities[key]['trait'][True] + p
        else:
            probabilities[key]['trait'][False] = probabilities[key]['trait'][False] + p
    
#    print("probabilities after update =",probabilities)
    
#    raise NotImplementedError


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    
    # for each 'gene' and 'trait' of a key (person) of probabilities, we'll normalize to 1
    
    for key in probabilities:
        normalizerGene = probabilities[key]['gene'][0] + probabilities[key]['gene'][1] + probabilities[key]['gene'][2]
        probabilities[key]['gene'][0] = probabilities[key]['gene'][0]/normalizerGene
        probabilities[key]['gene'][1] = probabilities[key]['gene'][1]/normalizerGene
        probabilities[key]['gene'][2] = probabilities[key]['gene'][2]/normalizerGene
    
        normalizerTrait = probabilities[key]['trait'][True] + probabilities[key]['trait'][False]
        probabilities[key]['trait'][True] = probabilities[key]['trait'][True]/normalizerTrait
        probabilities[key]['trait'][False] = probabilities[key]['trait'][False]/normalizerTrait        

    
    print("probabilities after normalization = ", probabilities)

if __name__ == "__main__":
    main()
