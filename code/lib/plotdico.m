function im = plotdico(dico)

    [D, N] = size(dico);

    n = ceil(sqrt(N));
    d = ceil(sqrt(D));

    im = 255 * ones(n*(d+1), n*(d+1));
    for k = 1:N
        i = mod(k - 1, n);
        j = floor((k - 1) / n);
        im((1:d) + i * (d+1), (1:d) + j * (d+1)) = reshape(dico(:,k),d,d);
    end
    im = im / 255;
    im(im < 0) = 0;
    im(im > 1) = 1;

    imshow(im);