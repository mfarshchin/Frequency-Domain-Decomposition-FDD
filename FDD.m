function [Frq,phi]=FDD(Input,Fs)

% Frequency Domain Decomposition (FDD) algorithm for modal analysis
% This code allows you to manually select the peaks by simply drawing a
% rectangle around the peaks.
% Programmer: Mohammad Farshchin, Ph.D candidate at The UofM
% Email: Mohammad.Farshchin@gmail.com
% Last modified: 4/6/2015

% Input: the name of input file that contains time history data
% Fs: sampling frequency
% Frq: identified frequencies
% phi: identified mode shapes
% Example: [Frq,phi]=FDD('Accelerations.xlsx',500);

% For detailed information about this method see: Brincker R, Zhang LM, Andersen P. Modal identification from ambient responses using Frequency Domain Decomposition. In: Proceedings of the 18th International Modal Analysis Conf., USA: San Antonio, 2000.

% -------------------------------------------------------------------------
% Initialization
close all
% -------------------------------------------------------------------------
% Import time history data: Processed accelereation data must be
% arranged in a columnwise format (one column for each measurement channel)
% Note that the acceleration data must be preprocessed (detrend, filtered etc.).
% Read acceleration data from the excel file
Acc=xlsread(Input);
display('FDD is in progress, please wait ...')
% -------------------------------------------------------------------------
% Compute Power Spectral Density (PSD) matrix.
% CPSD function, with default settings, is used to compute the cross power
% spectral density matrix. More sophisticated methods can also be
% applied for more accuracy.
for I=1:size(Acc,2)
    for J=1:size(Acc,2)
        [PSD(I,J,:),F(I,J,:)]=cpsd(Acc(:,I),Acc(:,J),[],[],[],Fs);
    end
end
Frequencies(:,1)=F(1,1,:);
% -------------------------------------------------------------------------
% Perform Modal Analysis (Use the Identifier function, below)
[Frq,phi,Fp,s1] = Identifier(PSD,Frequencies);
% Save results
save('IdResults.mat','phi','Fp','s1','Frequencies')
% -------------------------------------------------------------------------
% Print results
display('-----------------------------------------------------')
display('               Identification Results                ')
display('-----------------------------------------------------')
% Print frequencies
display('Identified frequencies')
for I=1:size(Frq,1)
    fprintf('Mode: %d; Modal Frequency: %6.4g (Hz)\n',I,Frq(I))
end
% Print Mode shapes
display('Related mode shapes')
for I=1:size(Frq,1)
    fprintf('Mode shape # %d:\n\n',I)
    disp(phi(:,I))
end

end

%% ---------------------------- subfunctions ------------------------------
function [Frq,phi,Fp,s1] = Identifier(PSD,F)

% Compute SVD of the PSD at each frequency
for I=1:size(PSD,3)
    [u,s,~] = svd(PSD(:,:,I));
    s1(I) = s(1);                                                          % First eigen values
    s2(I) = s(2,2);                                                        % Second eigen values
    ms(:,I)=u(:,1);                                                        % Mode shape
end

% Plot first singular values of the PSD matrix
figure
hold on
plot(F,mag2db(s1))
xlabel('Frequency (Hz)')
ylabel('1st Singular values of the PSD matrix (db)')

% -------------------------------------------------------------------------
% Peak selection
% a: Draw rectangles around peaks while holding left click
% b: Press "Space" key to continue peak selection
% c: Press "any other key" if you have selected a peak by mistake and want
% to ignore it
% -------------------------------------------------------------------------

Fp=[];% Frequencies related to selected peaks
NumPeaks=input('Enter the number of desired peaks:');
display('-----------------------------------------------------')
display('Peak selection procedure')
display('a: Draw rectangles around peaks while holding left click')
display('b: Press "Space" key to continue the peak selection')
display('c: Press "any other key" if you have selected a peak by mistake and want to ignore it')
k=0;
while k~=NumPeaks
    A=getrect;                                                                          % Draw a rectangle around the peak
    [~,P1]=min(abs(F-A(1)));
    [~,P2]=min(abs(F-(A(1)+A(3))));
    [~,B]=max(s1(P1:P2));
    Max=B+P1-1;                                                                         % Frequency at the selected peak
    scatter(F(Max),mag2db(s1(Max)),'MarkerEdgeColor','b','MarkerFaceColor','b')         % Mark this peak
    pause;key=get(gcf,'CurrentKey');
    Fp(end+1,:)=[Max,F(Max)];
    if strcmp(key,'space')
        % Press space to continue peak selection
        k=k+1;
        scatter(F(Max),mag2db(s1(Max)),'MarkerEdgeColor','g','MarkerFaceColor','g')      % Mark this peak as green
    else
        % Press any other key to ignore this peak
        Fp(end,:)=[];
        scatter(F(Max),mag2db(s1(Max)),'MarkerEdgeColor','r','MarkerFaceColor','r')      % Mark this peak as red
    end
end
% Number selected peaks, respectively
[~,Sr]=sort(Fp(:,2));
Fp=Fp(Sr,:);
clf
plot(F,mag2db(s1))
hold on
xlabel('Frequency (Hz)')
ylabel('1st Singular values of the PSD matrix (db)')
for I=1:size(Fp,1)
    scatter(Fp(I,2),mag2db(s1(Fp(I,1))),'MarkerEdgeColor','g','MarkerFaceColor','g')
    text(Fp(I,2), mag2db(s1(Fp(I,1)))*1.05, mat2str(I))
end
% Identified modal frequencies
Frq=Fp(:,2);
% Compute mode shapes for each selected peak
for J=1:size(Fp,1)
    [ug, ~, ~] = svd(PSD(:,:,Fp(J,1)));
    phi(:,J) = ug(:,1);
end

end






