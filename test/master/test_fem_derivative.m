A = rand(20); [dx dy dz] = derivative_3d(A); o = Test.d3(A,zeros(10,20),zeros(10,20),zeros(10,20),3); dz - o(3)

