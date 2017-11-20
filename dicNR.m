function dicNR()
	syms dx dy x0 y0 P P1 P2 P3 P4 P5 P6 P7 P8 P9 P10
	current_folder=pwd;
	addpath(strcat(current_folder,'\readimxstuff'));
	
	% save_as='Richard_CTC.mat';
	save_as='Richard_CTC_41_NR.mat';
	% image_count=max(size(FileName));

	inc=10;
	%define subset size
	subsize=41;
	stepsize=20;

	B=[(1+P2), P3, P1;
		P5, (P6+1), P4;
		0 0 1]
	X=[dx;dy;1];

	if exist(save_as,'file')
		load(save_as);
		% Proc.correlated_to=1;
		subpos=Proc.subpos;
		% guess_store=Proc.guess;
		guess_store=[0 0];
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		starting_subset=Proc.starting_subset;
		process_order=Proc.process_order;
		valid_subsets=Proc.valid_subsets;
		getIndex=Proc.getIndex;
		FileName=Proc.FileName;
		PathName=Proc.PathName;
		inc=Proc.inc;
		% Proc.inc=inc;
		X=Proc.WarpVec;
		B=Proc.Warp;
		stepsize=Proc.stepsize;
		subsize=Proc.subsize;
		current_image=Proc.correlated_to;
		% Proc.stepsize=stepsize;
		% Proc.subsize=subsize;
		% save_as='Richard_CTC_2_subset_41.mat';
		symbolic_warp(B,X)
	else 
		[FileName,PathName] = uigetfile('*.im7','Select the images','MultiSelect','on');
		Proc.FileName=FileName;
		Proc.PathName=PathName;
		Proc.inc=inc;
		Proc.stepsize=stepsize;
		Proc.subsize=subsize;
		Proc.Warp=B;
		Proc.WarpVec=X;
		symbolic_warp(B,X)
		for i=1:2
			% image_folder=strcat(PathName,FileName(i))
			image_folder = fullfile( PathName , FileName{i} );
			I{i}=readimx(image_folder);
		end
		F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
		mask=makeMask(F_in);
		% figure
		% imagesc(F_in)
		% polygon=impoly();
		% Proc.polygon=polygon;
		% mask=createMask(polygon);
		[subpos,mask_subsets,valid_subsets]=mask2subsets(mask,subsize,stepsize);
		Proc.subpos=subpos;
		Proc.mask_subsets=mask_subsets;
		Proc.valid_subsets=valid_subsets;
		[xguess,yguess,subx,suby]=seedPoints(PathName, FileName, subpos, mask_subsets,stepsize,subsize);
		Proc.guess=[xguess,yguess];
		guess_store=Proc.guess;
		guess=[guess_store(1),0,0,guess_store(2),0,0];
		starting_subset=whichSubpos(subpos,stepsize,subx,suby);
		Proc.starting_subset=starting_subset;
		[process_order,getIndex]=correlationOrder(subpos,starting_subset)
		Proc.process_order=process_order;
		Proc.getIndex=getIndex;
		Proc.correlated_to=1;
		current_image=1;
		save(save_as,'Proc')
	end

	image_folder = fullfile( PathName , FileName{1} );
	I{1}=readimx(image_folder);
	F_in=im2double(I{1}.Frames{1,1}.Components{1,1}.Planes{1,1});
    [r_F,c_F]=size(F_in);
 	image_count=max(size(FileName));

	
	% Proc.PathName=PathName;
	% Proc.FileName=FileName;
	% Proc.done=zeros([image_count,1]);
	% for i=1:image_count
	% 	fprintf('image %d\n',i);
	% 	if Proc.done(i)==0
	% 		image_folder = fullfile( PathName , FileName{i} );
	% 		G=readimx(image_folder);
	% 		tic
	% 		coef{i}=getBicubicValues(G.Frames{1,1}.Components{1,1}.Planes{1,1});
	% 		Proc.time(i)=toc;
	% 		Proc.coef{i}=coef{i};
	% 		Proc.done(i)=1;
	% 		save('Richard_CTC_coefficients.mat','Proc');
	% 	end
	% end


	elements=sum(sum(valid_subsets))
	size(process_order)

	for k=(current_image+inc):inc:image_count
		fprintf('image %d\n',k);
		tic
		Proc.im{k-1}.D=process_order;
		image_folder = fullfile( PathName , FileName{k} );
		I{3}=readimx(image_folder);
		G_in=im2double(I{3}.Frames{1,1}.Components{1,1}.Planes{1,1});
		if exist('coefficient_values.mat','file')
			load('coefficient_values.mat');
		else
			tic
			coef=getBicubicValues(G_in);
			toc
			save('coefficient_values.mat','coef');
			fprintf('calculated coefficients\n');
		end
		
		if (k==(1+inc))
			

			[PP(1,:),Corrr(1)]=NRtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(1,2),process_order(1,3)},'guess',guess,'coef',coef);
			for i=2:elements
				[PP(i,:),Corrr(i)]=NRtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',PP(i-1,:),'coef',coef);
			end
			for j=1:elements
				Proc.im{k-1}.D(j,6:11)=PP(j,:);
				Proc.im{k-1}.D(j,12)=Corrr(j);
			end
			
		else

			parfor i=1:elements
				[PP(i,:),Corrr(i)]=NRtracking2('undeformed image',F_in,'deformed image',G_in,'subset size',subsize,'subset position',subpos{process_order(i,2),process_order(i,3)},'guess',Proc.im{k-1-inc}.D(i,6:11),'coef',coef);
			end
			for j=1:elements
				Proc.im{k-1}.D(j,6:11)=PP(j,:);
				Proc.im{k-1}.D(j,12)=Corrr(j);
			end
		end

		toc
		Proc.correlated_to=k;
		save(save_as,'Proc');
	end


end