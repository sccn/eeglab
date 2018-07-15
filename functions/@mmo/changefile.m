% this function is called when the file is being saved

function obj = changefile(obj, newfile)

    movefile(obj.dataFile, newfile);
    obj.dataFile = newfile;
    obj.writable = false;
