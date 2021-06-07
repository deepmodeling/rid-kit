import re
import os
import shutil
import logging
import json
import datetime
import glob

iter_format = "%06d"
walker_format = "%03d"
task_format = "%02d"
log_iter_head = "iter " + iter_format + " task " + task_format + ": "


def make_iter_name(iter_index):
    return "iter." + (iter_format % iter_index)


def make_walker_name(walker_index):
    return (walker_format % walker_index)


def replace(file_name, pattern, subst):
    file_handel = open(file_name, 'r')
    file_string = file_handel.read()
    file_handel.close()
    file_string = (re.sub(pattern, subst, file_string))
    file_handel = open(file_name, 'w')
    file_handel.write(file_string)
    file_handel.close()


def create_path(path):
    if os.path.isdir(path):
        dirname = os.path.abspath(path)
        counter = 0
        while True:
            bk_dirname = dirname + ".bk%03d" % counter
            if not os.path.isdir(bk_dirname):
                os.system("mv {} {}".format(dirname, bk_dirname))
                print("{} has exists, a back-off {} generated!".format(path, bk_dirname))
                break
            counter += 1
    os.makedirs(path)


def copy_file_list(file_list, from_path, to_path):
    for jj in file_list:
        if os.path.isfile(os.path.join(from_path, jj)):
            shutil.copy(os.path.join(from_path, jj), to_path)
        elif os.path.isdir(os.path.join(from_path, jj)):
            cwd = os.getcwd()
            os.chdir(os.path.join(from_path, jj))
            files = glob.glob("*")
            os.chdir(cwd)
            os.makedirs(os.path.join(to_path, jj))
            for ff in files:
                shutil.copy(os.path.join(from_path, jj)+'/'+ff,
                            os.path.join(to_path, jj)+'/'+ff)


def checkfile(file_path):
    if os.path.exists(file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)


def log_task(message):
    header = repeat_to_length(" ", len(log_iter_head % (0, 0)))
    print(header + message)
    logging.info(header + message)


def log_iter(task, ii, jj):
    logging.info((log_iter_head + "%s") % (ii, jj, task))
    print((log_iter_head + "%s") % (ii, jj, task))


def repeat_to_length(string_to_expand, length):
    ret = ""
    for ii in range(length):
        ret += string_to_expand
    return ret


def cmd_append_log(cmd,
                   log_file):
    ret = cmd
    ret = ret + " 1> " + log_file
    ret = ret + " 2> " + log_file
    return ret


def print_list(tmp,
               suffix=""):
    mylist = ""
    for kk in tmp:
        if len(mylist) == 0:
            mylist = str(kk) + suffix
        else:
            mylist += "," + str(kk) + suffix
    # print(mylist)
    return mylist


def print_repeat_list(numb, item):
    mylist = ""
    for ii in range(numb):
        if ii == 0:
            mylist = str(item)
        else:
            mylist += "," + str(item)
    # print(mylist)
    return mylist


def record_iter(record, ii, jj):
    with open(record, "a") as frec:
        frec.write("%d %d\n" % (ii, jj))

# def record_task(record_file, iter_idx, task_idx):
#     if os.path.basename(record_file) == '':
#         print("please assign a valid record file path.")
#         raise RuntimeError

#     if os.path.exists(record_file):
#         with open(record_file, 'a') as record:
#             record.write("iteration: {} task: {}    task finished at {}\n".format(iter_idx, task_idx, datetime.datetime.now().strftime("%H:%M:%S %D")))
#     else:
#         with open(record_file, 'w') as record:
#             record.write("iteration: {} task: {}    task finished at {}\n".format(iter_idx, task_idx, datetime.datetime.now().strftime("%H:%M:%S %D")))
#     pass

# def get_checkpoint(record_file):
#     checkpoint = [-1, -1]
#     if not os.path.exists(record_file):
#         print("no record file exists, a new iteration will start.")
#         return [-1, -1]
#     else:
#         with open(record_file, 'r') as record:
#             for line in record.readlines():
#                 content = line.split()
#                 if len(content) == 0:
#                     continue
#                 else:
#                     checkpoint = [int(content[1]), int(content[3])]
#         print("The process will start after iteration {} task {}".format(checkpoint[0], checkpoint[1]))
#         return checkpoint


def get_checkpoint(record_file):
    checkpoint = [-1, -1]
    if not os.path.exists(record_file):
        print("no record file exists, a new iteration will start.")
        return [-1, -1]
    else:
        with open(record_file, 'r') as record:
            for line in record.readlines():
                content = line.split()
                if len(content) == 0:
                    continue
                else:
                    checkpoint = [int(content[0]), int(content[1])]
        print("The process will start after iteration {} task {}".format(
            checkpoint[0], checkpoint[1]))
        return checkpoint


if __name__ == '__main__':
    ret = print_list(['we', 'rng', 'edg'], suffix='0')
