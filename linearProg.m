function result = linearProg(blocks, rows, cols)
	
	blocksize = size(blocks, 1);
	n = size(blocks, 4);
	D = zeros(n, n, 4);		% compatibility of pieces (pairwise)

	B = 1e-6*[0, 0, 1; 1, 1, 1; 1, 0, 0];
	B = repmat(B, 1, 1, n);
	
	%%%%%%%%%%%%%%%%%%%%%%%% working on edges
	for o = 1:4
		blk = rot90(blocks, o-2);
		GiL = blk(:, end, :, :) - blk(:, end-1, :, :);	% an array of gradients for each channel
		GiL = reshape(GiL, blocksize, 3, n);
		GiL = [GiL; B];       
		muiL = mean(GiL);
		invSiL = zeros(3, 3, n);	% 3*3 covariance matrix
		for k = 1:n
			SiL = (GiL(:, :, k)-muiL(:, :, k)).'*(GiL(:, :, k)-muiL(:, :, k));
			invSiL(:, :, k) = inv(SiL);
		end
		for i = 1:n
			% the gradient from the right side of piece xi to the left side of piece xj
			GijLR = blk(:, 1, :, :) - blk(:, end, :, i);	
			GijLR = reshape(GijLR, blocksize, 3, n);
			for j = 1:n
				temp = GijLR(:, :, j) - muiL(:, :, i);
				D(i, j, o) = trace( temp * invSiL(:, :, i) * temp.' );
			end
		end
	end

	Dtemp = D;
	D(:, :, 1) = (Dtemp(:, :, 1) + Dtemp(:, :, 3).') / 2;
	D(:, :, 2) = (Dtemp(:, :, 2) + Dtemp(:, :, 4).') / 2;
	D(:, :, 3) = D(:, :, 1).';
	D(:, :, 4) = D(:, :, 2).';
	%%%%%%%%%%%%%%%%%%%%%%%% end of working of edges
	
	%%%%%%%%%%%%%%%% Converting MGC distance to weight
	n = size(D,  1);
	D = D + 1e-3*min(D(:))*rand(size(D));

	Dsortcol = sort(D);
	vminD = repmat(Dsortcol(1, :, :), n, 1, 1);
	vminD(vminD==D) = Dsortcol(2, :, :);

	Dtrans = permute(D, [2 1 3]);
	Dsortrow = sort(Dtrans);
	hminD = repmat(Dsortrow(1, :, :), n, 1, 1);
	hminD(hminD==Dtrans) = Dsortrow(2, :, :);
	hminD = permute(hminD, [2 1 3]);

	w = min(vminD, hminD) ./ D;
	%%%%%%%%%%%%%%%% End of MGC to weight
	
	%%%%%%%%%%%%%%%% Matching Constraints
	delta_x = [0; -1; 0; 1];
	delta_y = [1; 0; -1; 0];

	U = ones(n, n, 4);	% the universe of all possible oriented matches (i; j; o) between pieces

	[~, indj] = min(D.*U, [], 2);
	i = repmat((1:n)', 4, 1);
	j = indj(:);
	o = reshape(repmat(1:4, n, 1), [], 1);
	indA = sub2ind([n, n, 4], i, j, o);

	%%%%%%%%%%%%%%%% Minimizing the objective
	cvx_begin
		variables x(n) hx(4*n);
		minimize ( w(indA).' * hx );
		subject to
			-hx <= x(i) - x(j) - delta_x(o) <= hx;
			1 <= x <= cols;
	cvx_end

	cvx_begin
		variables y(n) hy(4*n);
		minimize ( w(indA).' * hy );
		subject to
			-hy <= y(i) - y(j) - delta_y(o) <= hy;
			1 <= y <= rows;
	cvx_end
	%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%% Constructing matrix A from matrix U
	subA = [i, j, o];
	subR = subA(max(hx, hy)>1e-5, :);
	indR = sub2ind([n, n, 4], subR);
	U(indR) = 0;
	%%%%%%%%%%%%%%%%

	%%%%%%%%%%%%%%%% Constructing the final image
	blocksize = size(blocks, 1);
	n = size(blocks, 4);

	x = x - min(x(:));
	x = round(x);
	y = y - min(y(:));
	y = round(y);

	result = ones((max(y)+1)*blocksize, (max(x)+1)*blocksize, 3);

	for k = 1:n
		rangecol = blocksize*x(k)+1 : blocksize*(x(k)+1);
		rangerow = blocksize*y(k)+1 : blocksize*(y(k)+1);
		result(rangerow, rangecol, :) = blocks(:, :, :, k);
		imshow(result);
	end
	title('The Final Result');
	%%%%%%%%%%%%%%%%
	
end
