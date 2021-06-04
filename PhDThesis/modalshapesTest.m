k = 1/44100;
c = 1470;
h = c*k;
N = 1/h;

Dxx = 1/h^2 * toeplitz([ -2, 1, zeros(1, N-2)]);

%% one step
Q = [2 * eye(N) + c^2 * k^2 * Dxx, -eye(N);
     eye(N), zeros(N)];

[phiQ, z] = eig(Q, 'vector');

s = log(z)/k;
s = imag(s);
[s, order] = sort(s);
phiQ = phiQ(:, order(s>0));

subplot(1,2,1); 
imagesc(real(phiQ)); 
title("Real part of $\phi$", 'interpreter', 'latex')
set(gca, 'Fontsize', 18)

subplot(1,2,2); 
imagesc(imag(phiQ)); 
title("Imaginary part of $\phi$", 'interpreter', 'latex')
set(gca, 'Fontsize', 16)
phiFinal = zeros(N);
set(gcf, 'color', 'white')

phiFinal = (sign(real(phiQ(1:30,:))) .* abs(phiQ(1:30,:)));

% imagVal = zeros(1, N);
% for i = 1:N
%     imagVal(i) = imag(phiQ(1:N,i))' * imag(phiQ(1:N,i)) ...
%         > imag(phiQ(N+1:end,i))' * imag(phiQ(N+1:end,i));
%     if ~imagVal(i)
%         phiFinal(:,i) = sign(real(phiQ(1:N, i))) .* abs(phiQ(1:N, i));
%     else
%         phiFinal(:,i) = sign(real(phiQ(N+1:end, i))) .* abs(phiQ(N+1:end, i));
%     end
% end
% 
% for i = 1:N
%     if sum(phiFinal(:,i)) < 0
%         phiFinal(:,i) =  phiFinal(:,i) * -1;
%     end
% end


imagesc(phiFinal)
%% two step
[phi, solut] = eig(Dxx, 'vector');
for i = 1:N
    if phi(1,i) < 0
        phi(:,i) =  phi(:,i) * -1;
    end
    if phiFinal(1,i) < 0
        phiFinal(:,i) =  phiFinal(:,i) * -1;
    end
end
phi = fliplr(phi);

for i = 1:N
    hold off;
    plot(phiFinal(:,i)*sqrt(2))
%     plot(abs(phiQ(1:30, i))/sin(pi/4))
    hold on;
%     plot(abs(phi(:,i)))
    plot(phi(:,i));
    diffSave(i) = sum(phi(:,i) - phiFinal(:,i)*sqrt(2));
    drawnow;
    pause(.05);
end

omega = 2/k * asin (c * k / 2 * sqrt(eig(-Dxx)));