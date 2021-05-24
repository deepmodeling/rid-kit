from context import rid
from rid import enhcMD, resMD, train, gen_rid
from rid.lib.utils import record_iter, get_checkpoint, make_iter_name
import json, os, glob
import argparse

def main(out_dir, mol_dir, rid_json, machine_json, cv_file, init_model, record_name="record.txt"):
    out_dir = os.path.abspath(out_dir)
    mol_dir = os.path.abspath(mol_dir)
    rid_json = os.path.abspath(rid_json)
    cv_file = os.path.abspath(cv_file)
    fp = open(rid_json, 'r')
    jdata = json.load(fp)
    fp.close()
    record_file = out_dir + record_name
    checkpoint = get_checkpoint(record_file)
    max_tasks = 10
    number_tasks = 8
    iter_numb = int(jdata['numb_iter'])
    prev_model = init_model

    if sum(checkpoint) < 0:
        print("prepare gen_rid")
        gen_rid.gen_rid (out_dir, mol_dir, rid_json)
    
    for iter_idx in range(iter_numb):
        if iter_idx > 0 :
            prev_model = glob.glob (out_dir + "/" + make_iter_name(iter_idx-1) + "/02.train/*pb")
        for tag in range(number_tasks):
            if iter_idx * max_tasks + tag <= checkpoint[0] * max_tasks + checkpoint[1]:
                # print(iter_idx, tag, 'pass')
                continue
            elif tag == 0:
                print("prepare gen_enhc")
                enhcMD.make_enhc(iter_idx, rid_json, prev_model, mol_dir, cv_file ,base_dir=out_dir)
                record_iter(record_file, iter_idx, tag)
            elif tag == 1:
                print("prepare run enhc")
                enhcMD.run_enhc(iter_idx, rid_json, machine_json, base_dir=out_dir)
                record_iter(record_file, iter_idx, tag)
            elif tag == 2:
                print("prepare post enhc")
                enhcMD.post_enhc(iter_idx, rid_json, machine_json, base_dir=out_dir)
                record_iter(record_file, iter_idx, tag)
            elif tag == 3:
                print("prepare gen_res")
                resMD.make_res(iter_index=iter_idx, json_file=rid_json, cv_file=cv_file, mol_path=mol_dir, base_dir=out_dir)
                record_iter(record_file, iter_idx, tag)
            elif tag == 4:
                print("prepare run res")
                resMD.run_res (iter_index=iter_idx, json_file=rid_json, machine_json=machine_json, base_dir=out_dir)
                record_iter(record_file, iter_idx, tag)
            elif tag == 5:
                print("prepare post res")
                resMD.post_res(iter_index=iter_idx, json_file=rid_json, machine_json=machine_json, cv_file=cv_file, base_dir=out_dir)
                record_iter(record_file, iter_idx, tag)
            elif tag == 6:
                print("prepare gen train")
                train.make_train(iter_index=iter_idx, json_file=rid_json, base_dir=out_dir)
                record_iter(record_file, iter_idx, tag)
            elif tag == 7:
                print("prepare run train")
                train.run_train(iter_index=iter_idx, json_file=rid_json, machine_json=machine_json, cv_file=cv_file, base_dir=out_dir)
                record_iter(record_file, iter_idx, tag)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("JSON", type=str, 
                        help="The json parameter")
    parser.add_argument("-c", "--cv", type=str, 
                        help="The selected CV json.")
    parser.add_argument("-s", "--machine", type=str, 
                        help="The machine settings")        
    parser.add_argument("-i", "--mol", type=str, 
                        help="The initial conformation folder.")
    parser.add_argument("-o", "--out", type=str, 
                        help="The out path.")
    parser.add_argument("-m", "--models", default=[], nargs = '*', type=str, 
                        help="The init guessed model")    
    args = parser.parse_args()

    # print(args.JSON)
    # print(args.machine)
    # print(args.mol)
    # print(args.out)
    # print(args.models)
    # print(args.cv)

    main(out_dir=args.out, 
         mol_dir=args.mol, 
         rid_json=args.JSON, 
         machine_json=args.machine, 
         cv_file=args.cv, 
         init_model=args.models, 
         record_name="record.txt")
