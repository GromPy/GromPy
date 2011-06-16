from grompy import matrix,rvec
from grompy.types import XX,YY,ZZ

# from vec,h
def copy_rvec(a = rvec):
    b     = rvec()
    b[XX] = a[XX]
    b[YY] = a[YY]
    b[ZZ] = a[ZZ]

    return b
#static inline void copy_rvec(const rvec a,rvec b)
#{
#  b[XX]=a[XX];
#  b[YY]=a[YY];
#  b[ZZ]=a[ZZ];
#}

def copy_mat(a = matrix):
    b     = matrix()
    b[XX] = copy_rvec(a[XX])
    b[YY] = copy_rvec(a[YY])
    b[ZZ] = copy_rvec(a[ZZ])

    return b
#static inline void copy_mat(matrix a,matrix b)
#{
#  copy_rvec(a[XX],b[XX]);
#  copy_rvec(a[YY],b[YY]);
#  copy_rvec(a[ZZ],b[ZZ]);
#}