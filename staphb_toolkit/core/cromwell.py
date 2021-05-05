import urllib.request
import glob
import os,sys


def get_cromwell_engine(version_tag=""):
    lib_path = os.path.abspath(os.path.dirname(__file__) + '../' + '../' + '/lib')

    #download latest
    if not version_tag:
        #get url to latest version
        cromwell_latest_url = "https://github.com/broadinstitute/cromwell/releases/latest"
        version_tag = urllib.request.urlopen(cromwell_latest_url).geturl().split('/')[-1]

    cromwell_download_url = f"https://github.com/broadinstitute/cromwell/releases/download/{version_tag}/cromwell-{version_tag}.jar"
    download_path = os.path.abspath(os.path.join(lib_path,f"cromwell-{version_tag}.jar"))
    #download cromwell
    if not os.path.isfile(download_path):
        print("Updating to cromwell version:",version_tag)
        #remove old versions
        old_list = glob.glob(os.path.abspath(os.path.join(lib_path,"cromwell-*.jar")))
        for file in old_list:
            try:
                os.remove(file)
            except:
                print("Couldn't delete old cromwell version: ",file)

        #download new version
        with urllib.request.urlopen(cromwell_download_url) as response, open(download_path,'wb') as outfile:
            data = response.read()
            outfile.write(data)

        print("Cromwell sucessfully updated to:",version_tag)
