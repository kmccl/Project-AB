%creating a common reference; reference becomes last electrode
%function q = comreff(d)


function q = comreff(d)
[i j]=size(d);

if (size(d,2)==j)
	q = zeros(size(d,1),j+1);         
	r = mean(d');
	r = r';
	q = d - (r*ones(1,j));
	q(:,j+1) = -r;
end



