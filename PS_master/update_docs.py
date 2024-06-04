'''
Script to update the documentation website
'''

import os
import shutil
import filecmp
import yaml
import copy
import fileinput

# Defining machine
machine = 'PS'

# Ignore folders
folder_ignore_list = ['source']

# Extension of figures to include in the docs
delete_old_files = True
copy_images = True
fig_extension = ['png', 'jpg', 'html', 'PNG', 'JPG', 'HTML']

# Creating folders and comparing files already uploaded to the docs
base_imp_path = './impedance/'
base_doc_path = './mkdocs/'

if machine not in os.listdir(os.path.join(base_doc_path, 'docs')):
    if not os.path.exists(os.path.join(base_doc_path, 'docs', machine)):
        os.makedirs(os.path.join(base_doc_path, 'docs', machine))
else:
    if delete_old_files:
        for folder in os.listdir(os.path.join(base_doc_path, 'docs', machine)):
            shutil.rmtree(os.path.join(base_doc_path, 'docs', machine, folder),
                          ignore_errors=True)

# Loading existing nav file
with open(base_doc_path + '/mkdocs.yml') as nav_file:
    loaded_nav = yaml.safe_load(nav_file)
for index, machine_dict in enumerate(loaded_nav['nav']):
    if list(machine_dict.keys())[0] == machine:
        machine_dict = machine_dict
        machine_index = index
        break

new_nav = list(loaded_nav['nav'])


# Define capitalize function
def capwords(string):
    return ' '.join(['%s%s' % (element[0].upper(), element[1:]) for element in string.split(' ')])


# Treating main index file
path = os.path.join(base_imp_path, '..')
for name in os.listdir(path):
    if name.split('.')[-1] == 'md':
        main_md = name
        save_path = os.path.join(base_doc_path, 'docs', machine)
        nav_path = os.path.join(machine, name)
        if not os.path.exists(os.path.join(save_path, name)):
            shutil.copyfile(os.path.join(path, name),
                            os.path.join(save_path, name))
        else:
            if not filecmp.cmp(os.path.join(path, name),
                               os.path.join(save_path, name)):
                shutil.copyfile(os.path.join(path, name),
                                os.path.join(save_path, name))
        break


new_nav[machine_index][machine] = [{'Overview': nav_path}]
group_list = ['Overview']
eq_list = [[]]
sub_list = [[[]]]

# Looping over all folders and copying .md files

for path, subdirs, files in os.walk(base_imp_path):
    for name in files:

        if name.split('.')[-1] == 'md':
            path_split = path[len(base_imp_path):].split('/')

            group = path_split[0]
            group_renamed = capwords(group.replace('_', ' '))
            if group in folder_ignore_list:
                continue
            save_path = os.path.join(base_doc_path, 'docs', machine, group)
            if group_renamed not in group_list:
                group_list.append(group_renamed)
                index_group = len(group_list) - 1
                new_nav[machine_index][machine].append({group_renamed: []})
                eq_list.append([])
                sub_list.append([[]])
            else:
                for index, group_find in enumerate(group_list):
                    if group_renamed == group_find:
                        index_group = index
                        break
            nav_path = os.path.join(machine, group, name)
            if not os.path.exists(save_path):
                os.makedirs(save_path)

            if len(path_split) > 1:

                equipment = path_split[1]
                eq_renamed = capwords(equipment.replace('_', ' '))
                if equipment in folder_ignore_list:
                    continue
                save_path = os.path.join(base_doc_path, 'docs', machine,
                                         group, equipment)
                if eq_renamed not in eq_list[index_group]:
                    eq_list[index_group].append(eq_renamed)
                    index_eq = len(eq_list[index_group]) - 1
                    new_nav[machine_index][machine][index_group][group_renamed].append({
                                                                                       eq_renamed: []})
                    sub_list[index_group].append([])
                else:
                    for index, eq_find in enumerate(eq_list[index_group]):
                        if eq_renamed == eq_find:
                            index_eq = index
                            break
                nav_path = os.path.join(machine, group, equipment, name)
                if not os.path.exists(save_path):
                    os.makedirs(save_path)
            else:
                pass

            if len(path_split) > 2:

                subcase = path_split[2]
                sub_renamed = capwords(subcase.replace('_', ' '))
                if subcase in folder_ignore_list:
                    continue
                save_path = os.path.join(base_doc_path, 'docs', machine,
                                         group, equipment, subcase)
                if sub_renamed not in sub_list[index_group][index_eq]:
                    sub_list[index_group][index_eq].append(sub_renamed)
                    index_sub = len(sub_list[index_group][index_eq]) - 1
                    new_nav[machine_index][machine][index_group][group_renamed][index_eq][eq_renamed].append({
                                                                                                             sub_renamed: []})
                else:
                    for index, sub_find in enumerate(sub_list[index_group][index_eq]):
                        if sub_renamed == sub_find:
                            index_sub = index
                            break
                nav_path = os.path.join(machine, group, equipment,
                                        subcase, name)
                if not os.path.exists(save_path):
                    os.makedirs(save_path)
            else:
                sub_renamed = ''

            with fileinput.FileInput(os.path.join(path, name), inplace=True) as file:
                math_bool = False
                for line in file:
                    if '```math' in line:
                        math_bool = True
                    new_line = line.replace('```math', '$$')
                    new_line = new_line.replace('$`', '$')
                    new_line = new_line.replace('`$', '$')
                    if math_bool and ('```' in new_line):
                        new_line = new_line.replace('```', '$$')
                        math_bool = False
                    print(new_line, end='')

            if not os.path.exists(os.path.join(save_path, name)):
                shutil.copyfile(os.path.join(path, name),
                                os.path.join(save_path, name))
            else:
                if not filecmp.cmp(os.path.join(path, name),
                                   os.path.join(save_path, name)):
                    shutil.copyfile(os.path.join(path, name),
                                    os.path.join(save_path, name))

            if len(path_split) == 1:
                nav_path = os.path.join(machine, group, name)
                new_nav[machine_index][machine][index_group][group_renamed].append(
                    {'Overview': nav_path})
                eq_list[index_group].append([])
            elif len(path_split) == 2:
                nav_path = os.path.join(machine, group, equipment, name)
                new_nav[machine_index][machine][index_group][group_renamed][index_eq][eq_renamed].append({
                                                                                                         'Overview': nav_path})
                sub_list[index_group][index_eq].append([])
            elif len(path_split) == 3:
                nav_path = os.path.join(
                    machine, group, equipment, subcase, name)
                new_nav[machine_index][machine][index_group][group_renamed][index_eq][eq_renamed][index_sub][sub_renamed] = nav_path

# Sorting
sorted_nav = copy.deepcopy(new_nav)
group_list_sorted = copy.deepcopy(group_list)
eq_list_sorted = copy.deepcopy(eq_list)
for index_group in range(1, len(group_list)):
    sorted_index_gr = sorted(group_list[1:]).index(
        group_list[1:][index_group - 1]) + 1
    sorted_nav[machine_index][machine][sorted_index_gr] = copy.deepcopy(
        new_nav[machine_index][machine][index_group])
    group = group_list[index_group]
    group_list_sorted[sorted_index_gr] = copy.deepcopy(group)
    eq_list_sorted[sorted_index_gr] = copy.deepcopy(eq_list[index_group])
    if (len(eq_list[index_group]) > 0) and (eq_list[index_group][0] == []):
        start_point = 1
    else:
        start_point = 0
    for index_eq in range(start_point, len(eq_list[index_group][start_point:]) + start_point):
        sorted_index_eq = sorted(eq_list[index_group][start_point:]).index(
            eq_list[index_group][start_point:][index_eq - start_point]) + start_point
        sorted_nav[machine_index][machine][sorted_index_gr][group][sorted_index_eq] = copy.deepcopy(
            new_nav[machine_index][machine][index_group][group][index_eq])
        equipment = eq_list[index_group][index_eq]
        eq_list_sorted[sorted_index_gr][sorted_index_eq] = copy.deepcopy(
            equipment)

# Reducing nav size
reduced_nav = copy.deepcopy(sorted_nav)
for index_group in range(1, len(reduced_nav[machine_index][machine])):
    if len(reduced_nav[machine_index][machine][index_group][group_list_sorted[index_group]]) == 1:
        reduced_nav[machine_index][machine][index_group][group_list_sorted[index_group]
                                                         ] = reduced_nav[machine_index][machine][index_group][group_list_sorted[index_group]][0]['Overview']
    else:
        for index_eq in range(len(reduced_nav[machine_index][machine][index_group][group_list_sorted[index_group]])):
            if 'Overview' in reduced_nav[machine_index][machine][index_group][group_list_sorted[index_group]][index_eq].keys():
                continue
            if len(reduced_nav[machine_index][machine][index_group][group_list_sorted[index_group]][index_eq][eq_list_sorted[index_group][index_eq]]) == 1:
                reduced_nav[machine_index][machine][index_group][group_list_sorted[index_group]][index_eq][eq_list_sorted[index_group][index_eq]] = \
                    reduced_nav[machine_index][machine][index_group][group_list_sorted[index_group]
                                                                     ][index_eq][eq_list_sorted[index_group][index_eq]][0]['Overview']

# Writing in nav file
with open(base_doc_path + '/mkdocs.yml', 'w') as nav_file:
    final_nav = dict(loaded_nav)
    final_nav['nav'] = reduced_nav
    yaml.dump(final_nav, nav_file, sort_keys=False)


# Looping over folders and copying images
if copy_images:
    for path, subdirs, files in os.walk(base_imp_path):
        for name in files:
            if name.split('.')[-1] in fig_extension:
                path_split = path[len(base_imp_path):].split('/')
                save_path = os.path.join(
                    base_doc_path, 'docs', machine, *path_split)
                if not os.path.exists(save_path):
                    os.makedirs(save_path)
                if not os.path.exists(os.path.join(save_path, name)):
                    shutil.copyfile(os.path.join(path, name),
                                    os.path.join(save_path, name))
                else:
                    if not filecmp.cmp(os.path.join(path, name),
                                       os.path.join(save_path, name)):
                        shutil.copyfile(os.path.join(path, name),
                                        os.path.join(save_path, name))
