IDX_ok_1 = find(allvarf<1.0);
IDX_ok_2 = find(allvarf<2.0);
IDX_ok_3 = find(allvarf<3.0);
IDX_ok_4 = find(allvarf<4.0);
allESTPS_1= allESTPS(:,IDX_ok_1);
allESTPS_2= allESTPS(:,IDX_ok_2);
allESTPS_3= allESTPS(:,IDX_ok_3);
allESTPS_4= allESTPS(:,IDX_ok_4);


def_min_1=min(allESTPS(2,IDX_ok_1));
def_max_1=max(allESTPS(2,IDX_ok_1));
def_mean_1=mean(allESTPS(2,IDX_ok_1));
def_std_1=std(allESTPS(2,IDX_ok_1));

dem_min_1=min(allESTPS(1,IDX_ok_1));
dem_max_1=max(allESTPS(1,IDX_ok_1));
dem_mean_1=mean(allESTPS(1,IDX_ok_1));
dem_std_1=std(allESTPS(1,IDX_ok_1));

def_min_2=min(allESTPS(2,IDX_ok_2));
def_max_2=max(allESTPS(2,IDX_ok_2));
def_mean_2=mean(allESTPS(2,IDX_ok_2));
def_std_2=std(allESTPS(2,IDX_ok_2));


dem_min_2=min(allESTPS(1,IDX_ok_2));
dem_max_2=max(allESTPS(1,IDX_ok_2));
dem_mean_2=mean(allESTPS(1,IDX_ok_2));
dem_std_2=std(allESTPS(1,IDX_ok_2));

def_min_3=min(allESTPS(2,IDX_ok_3));
def_max_3=max(allESTPS(2,IDX_ok_3));
def_mean_3=mean(allESTPS(2,IDX_ok_3));
def_std_3=std(allESTPS(2,IDX_ok_3));


dem_min_3=min(allESTPS(1,IDX_ok_3));
dem_max_3=max(allESTPS(1,IDX_ok_3));
dem_mean_3=mean(allESTPS(1,IDX_ok_3));
dem_std_3=std(allESTPS(1,IDX_ok_3));

def_min_4=min(allESTPS(2,IDX_ok_4));
def_max_4=max(allESTPS(2,IDX_ok_4));
def_mean_4=mean(allESTPS(2,IDX_ok_4));
def_std_4=std(allESTPS(2,IDX_ok_4));


dem_min_4=min(allESTPS(1,IDX_ok_4));
dem_max_4=max(allESTPS(1,IDX_ok_4));
dem_mean_4=mean(allESTPS(1,IDX_ok_4));
dem_std_4=std(allESTPS(1,IDX_ok_4));