#!/usr/bin/env python3

import os,re

#replace absolute and relative paths with paths in docker container
def path_replacer(args,cwd):
    #generate mapping of paths and container paths
    path_map = {cwd:'/data'} #{"/path/outside":"/path/incontainer"}
    if not args:
        return "",path_map
    #search pattern for absolute and relative paths
    regex_path_pattern = '^[\/,~\/,.]\S*'
    counter = 1
    arg_string = '' #string of final arguments for docker command
    for arg in args:
        #check if argument is a path
        re_obj = re.match(regex_path_pattern,arg)
        #if it is add a path mapping
        if re_obj:
            path = re_obj.group()
            basename = os.path.basename(path)
            dirname = os.path.dirname(path)
            #check if we have already created a mount point for this location
            if dirname in path_map.keys():
                mountname = os.path.dirname(path_map[dirname])
                path_map[dirname] = mountname+'/'
            else:
                path_map[dirname] = '/mount'+str(counter)+'/'
                counter += 1
            arg_string = arg_string + ' ' + path_map[dirname]+'/'+basename
        #if it's not add the argument to final string
        else:
            arg_string = arg_string + ' ' + arg

    return arg_string,path_map
