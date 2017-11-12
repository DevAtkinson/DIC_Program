function NR(subsize,numP,A,X)
	F = sym('F%d%d', [subsize subsize]);
	G = sym('G%d%d', [subsize subsize]);
	P = sym('P%d', [numP]);
	syms dx dy;
	W=A*X;
	[r_g,c_g]=size(G);
	[Xmesh,Ymesh]=meshgrid(1:1:c_g,1:1:r_g);
	% Int = (Xmesh',Ymesh',G','cubic')

end