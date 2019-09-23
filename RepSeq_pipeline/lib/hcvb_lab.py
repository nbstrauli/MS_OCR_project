def is_bcr(file_name):
    uc_file_name = file_name.upper()
    indicators = ['BCR',  'IGM', 'IGG', 'IGH', 'IGVH']
    #TODO: add key work for TCR and check the key work is not in filename
    bcr = any([ x in uc_file_name for x in indicators])
    return bcr

def raw_data_dir_to_analysi_dir(raw_data_dir):
    """ remove the 'Auto_user_SN2 in dirname
    """
    return raw_data_dir.lstrip('Auto_user_SN2-')

def analysis_dir_to_raw_data_dir(analysis_dir):
    """we striped the Auto_user part in raw_data_dir_to_analysi_dir
    """
    raw_data_dir= os.path.basename(analysis_dir)
    if os.path.basename(raw_data_dir)[0].isdigit():
        raw_data_dir =  'Auto_user_SN2-' + raw_data_dir
    return raw_data_dir

def extract_bcr_id(item):
    csv= item.split('/')[-1]
    if csv.endswith('.csv') and is_bcr(csv):
        sample_name = csv.rstrip('.csv')
        return '_'.join(sample_name.split('_')[0:2])

def get_changed_bcr_ids(copy_log_of_today):
     ids = []
     with open(copy_log_of_today, 'r') as f:
         for item in f:
             current_id = extract_bcr_id(item.rstrip()[24:])
             if current_id:
                ids.append(current_id)
                ids.append(current_id.split('_')[0])
     return(list(set(ids)))

def extract_tcr_id(item):
    csv= item.split('/')[-1]
    if csv.endswith('.csv') and not is_bcr(csv):
        sample_name = csv.rstrip('.csv')
        return '_'.join(sample_name.split('_')[0:2])

def get_changed_tcr_ids(copy_log_of_today):
     ids = []
     with open(copy_log_of_today, 'r') as f:
         for item in f:
             current_id = extract_tcr_id(item.rstrip()[24:])
             if current_id:
                ids.append(current_id)
                ids.append(current_id.split('_')[0])
     return(list(set(ids)))


def get_changed_run(copy_log_of_today):
     with open(copy_log_of_today, 'r') as f:
         for item in f:
             if 'copying' in item:
                 run_path = os.path.dirname(item.split()[-1])
                 return os.path.dirname(item.split()[-1]).split('/')[4]
