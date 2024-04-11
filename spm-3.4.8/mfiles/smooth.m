function B = smooth(A)

[ny,nx]=size(A);

B = A;
B(2:ny-1,2:nx-1) = B(2:ny-1,2:nx-1)+3*A(2:ny-1,2:nx-1)+A(2:ny-1,1:nx-2)+A(3:ny,1:nx-2)+A(3:ny,2:nx-1)+A(3:ny,3:nx)+A(2:ny-1,3:nx)+A(1:ny-2,3:nx)+A(1:ny-2,2:nx-1)+A(1:ny-2,1:nx-2);
B(2:ny-1,2:nx-1) = B(2:ny-1,2:nx-1)/12;