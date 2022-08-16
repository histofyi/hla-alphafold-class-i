from common import perform_action
from functions import generate_rmsd_set


errors = []

listname = 'human_class_i.csv'

errors = perform_action(generate_rmsd_set, listname)

print (errors)




























