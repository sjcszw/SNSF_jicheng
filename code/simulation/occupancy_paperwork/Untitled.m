
A1 = model.A
B1 = model.B
C1 = model.C
D1 = model.D

i = 1;
x1 = A1* x_cl(:,i) + B1*u_cl(:,i);
y1 = C1*x1 + D1*u_cl(:,i+1);
x1- x_cl(:,i+1)
y1 - y_cl(:,i+1)

x2 = A* x_cl(:,i) + B*u_cl(:,i);
y2 = C*x2 + D*u_cl(:,i+1);
x2- x_cl(:,i+1)
y2 - y_cl(:,i+1)


C*x_cl(:,i+1) + D*u_cl(:,i+1) - y_cl(:,i+1)