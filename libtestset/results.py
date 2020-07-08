import csv


def write_results(options, systems):
    pass

# Minimal example for later
def proof_of_concept():
    categories = ["REACTION", "REFERENCE", "DFTB3"]
    data = [
        {"REACTION": "2 A + B -> 2 C", "REFERENCE": 103.910, "DFTB3": 111.002},
        {"REACTION": "3 A + B + C -> D", "REFERENCE": 150.070, "DFTB3": 155.801}
    ]
    csv_file = "output.csv"

    with open(csv_file, "w") as out:
        writer = csv.DictWriter(out, fieldnames=categories)
        writer.writeheader()
        for entry in data:
            writer.writerow(entry)
