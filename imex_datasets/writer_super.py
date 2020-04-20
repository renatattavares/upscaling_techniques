with open('original_files/permi_final.txt') as permi:
    permx = permi.readlines()

with open('original_files/permj_final.txt') as permj:
    permy = permj.readlines()

with open('original_files/permk_final.txt') as permk:
    permz = permk.readlines()

with open('super.dat') as super:
    dataset = super.readlines()

i = 0

for line in dataset:
    if '*PERMI *ALL' in line:
            print('Here {}'.format(i))

    if '*MODEL *OILWATER' in line:
            print('Here {}'.format(i))

    i += 1

heading = dataset[:224486]
final = dataset[897692:]

with open('new_super.dat', 'w') as new:
    new.writelines(heading)
    new.write('\n*PERMI *ALL     ** An array of geostatistically distributed\n')
    new.writelines(permx)
    new.write('\n')
    new.write('\n*PERMJ *ALL     ** An array of geostatistically distributed\n')
    new.writelines(permy)
    new.writelines('\n')
    new.write('\n*PERMK *ALL     ** An array of geostatistically distributed\n')
    new.writelines(permz)
    new.write('\n')
    new.writelines(final)
