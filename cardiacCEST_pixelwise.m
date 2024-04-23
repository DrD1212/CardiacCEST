%For 2D CEST data (acquired via interleave mode)

%Calculation of CEST spectrum pixel by pixel including B0 correction to create asymmetry maps & to quantify the CEST effect
%Requires brukertools, esp. "readBrukerParamFile" and "readBruker2dseq" to load 2dseq MRI data

%Possibility to compare pre and post infusion data
%Possibility to filter data via DOSE method

%The file path and the quantification parameters need to be adjusted manually

%Spectrum for each pixel can then be found in "Conclusion_frequencies" or "Conclusion_data".
%Asymmetry values at a certain offset can be found in "Conclusion_asym".

%Main part: Calculation of pixel-by-pixel CEST and associated asymmetry spectra (with DOSE filtering if necessary)

%Add 1: Calculation of a global CEST spectrum (unfiltered)
%Add 2: Calculation of pixel-by-pixel CEST maps (unfiltered)
%Add 3: Calculation of a global asymmetry spectrum (unfiltered)
%Add 4: Processing of DOSE-filtered CEST data
%Add 5: Exclusion of outliers in DOSE-filtered data & calculation of filtered global CEST and asymmetry spectrum


%%
%clear all
%close all
clearvars -except ROI copy_ROI_Asymspectrum A_Prae_Ergebnisse copy_ROI_Asymspectrum_korr copy_ROI_Zspectrum_korr copy_ROI_Zspectrum_freq_korr ababa Mittelungsspeicher Anzahlspeicher dataset

%%
%Load data via bruker function (adjust path for CEST and S0 data!)

%CEST data:
visu_CEST = readBrukerParamFile('C:\Users\d_scha31\Desktop\CEST\in vivo Messdaten\Herz-CEST-Ajay\11252_20230512_D1\48\pdata\1\visu_pars');
scan = readBruker2dseq('C:\Users\d_scha31\Desktop\CEST\in vivo Messdaten\Herz-CEST-Ajay\11252_20230512_D1\48\pdata\1\2dseq',visu_CEST);

%S0 data:
visu_S0 = readBrukerParamFile('C:\Users\d_scha31\Desktop\CEST\in vivo Messdaten\Herz-CEST-Ajay\11252_20230512_D1\47\pdata\1\visu_pars');
scanS0 = readBruker2dseq('C:\Users\d_scha31\Desktop\CEST\in vivo Messdaten\Herz-CEST-Ajay\11252_20230512_D1\47\pdata\1\2dseq',visu_S0);


AX=squeeze(scan);
BX=squeeze(scanS0);

CESTdata_3D = permute(AX, [3 1 2]); 
S0data_3D = permute(BX, [3 1 2]);

A=size(CESTdata_3D); 
B=size(S0data_3D);


%%
%%Adjust parameters:
Range=2000; 
Int_Offset=512; %Single offset for quantification
lowBound_AUC=200; %AuC for quantification
highBound_AUC=800;


PraeMeas=0; %0 = pre; 1= post
ROI_exist=0; %0 = draw new ROI; 1 = maintain ROI
Seperate=1; %0 = without seperation of data; 1 = seperation of data (required for DOSE filtering)


MinCluster=1; %Cluster size

%%
%%Definition of parameters / matrices:
NROffsets=A(1); 
NROffsets_S0=B(1); 
Pixelx=A(2);
Pixely=A(3); 


frequenzen=linspace(-Range,Range,NROffsets); 

Conclusion_frequencies=zeros(Pixelx,Pixely,NROffsets); 
Conclusion_data=zeros(Pixelx,Pixely,NROffsets); 
Conclusion_asym=zeros(Pixelx,Pixely,1);

Conclusion_asymKurven=zeros(Pixelx,Pixely,Range);
Conclusion_asymfrequencies=zeros(Pixelx,Pixely,Range);

Conclusion_AUC=zeros(Pixelx,Pixely,1);

Conclusion_Minimum=zeros(Pixelx,Pixely,1);

Conclusion_spline_x=zeros(Pixelx,Pixely,(2*Range+1));
Conclusion_spline_y=zeros(Pixelx,Pixely,(2*Range+1));


Conclusion_seperate_x=zeros(Pixelx,Pixely,Range);
Conclusion_seperate_y=zeros(Pixelx,Pixely,Range);
Conclusion_Center=zeros(Pixelx,Pixely,1);

Conclusion_sat=zeros(Pixelx,Pixely,1);


CESTdata_rot=squeeze(CESTdata_3D(:,:,:));
S0data_rot=squeeze(S0data_3D(:,:,:));

CESTdata= permute(CESTdata_rot,[1 3 2]); 
S0data= permute(S0data_rot,[1 3 2]);


%%
%Pre-calculation:

for i=1:NROffsets
   data(i,:,:)= squeeze(CESTdata(i,:,:)); 
end

for i=1:NROffsets_S0
   data_S0(i,:,:)= squeeze(S0data(i,:,:)); 
end


%Insert: Frequency list (as data are acquired in interleave mode)----------
steps=NROffsets;
stepsize=Range/((steps-1)/2);

Offsets=zeros(steps,1);

%left side
for i=1:((steps-1)/2)
   Offsets((2*i-1),1)=-Range+(i-1)*stepsize; 
end

%right side
for i=1:((steps-1)/2)
   Offsets((2*i),1)=Range-(i-1)*stepsize; 
end

Offsets(steps,1)=0;

%End of the insert---------------------------------------------------------


%Define ROI
AA=squeeze(CESTdata(1,:,:));
if(ROI_exist==0)
figure();
imagesc(AA(:,:));
ROI=roipoly;
end

%%

%Pixelwise CEST calculation (for loop; calculation for all pixels within the ROI)

for  k=1:Pixely
for  j=1:Pixelx

number = sum(ROI(:)~=0); 

if(ROI(j,k)==1)
clearvars -except Conclusion_frequencies Conclusion_data Conclusion_asym Conclusion_asymKurven Conclusion_asymfrequencies Conclusion_Minimum Conclusion_AUC Int_Offset CESTpath CESTdatapath CESTdata handles.CEST.CEST NROffsets NROffsets_S0 Pixelx Pixely Range frequenzen data data_S0 Offsets k j asym_Bound MinCluster ROI_exist PraeMeas ROI highBound_AUC lowBound_AUC copy_ROI_Asymspectrum A_Prae_Ergebnisse copy_ROI_Asymspectrum_korr copy_ROI_Zspectrum_korr copy_ROI_Zspectrum_freq_korr AA ababa Conclusion_spline_x Conclusion_spline_y Conclusion_seperate_x Conclusion_seperate_y Conclusion_Center Seperate Mittelungsspeicher Anzahlspeicher dataset Seperate_quant B Conclusion_sat
if (NROffsets_S0 > 1)
for n=1:NROffsets 
    spectrum_1(n)=data(n,j,k)/((data_S0(1,j,k)+data_S0(2,j,k))/2); 
end
end

if (NROffsets_S0 == 1)  
for n=1:NROffsets
    spectrum_1(n)=data(n,j,k)/((data_S0(1,j,k))); 
end
end

%Sort data (as data are acquired in interleave mode)
clear i
for i=1:((NROffsets-1)/2) %vorher 31
   CESTdata_2(1,i)=spectrum_1(1,(2*i-1)); 
end

for i=1:((NROffsets-1)/2)
    CESTdata_2(1,i+((NROffsets+1)/2))=spectrum_1(1,(NROffsets+1-2*i)); 
end

CESTdata_2(1,((NROffsets+1)/2))=spectrum_1(1,NROffsets); 


%%

%B0-correction (based on data spline)
freqRange=linspace(-Range,Range,size(CESTdata_2,2)); 

freqRangeHR=linspace(-Range,Range,(2*Range)+1); 

pros = spline(freqRange(1:end),CESTdata_2(1:end),freqRangeHR);

%Find the minimum of the spline
[min_val,idx]=min(pros(:)); 
[row,col]=ind2sub(size(pros),idx);
Minimum=((2*Range/(2*Range))*col)-Range; 
fprintf('\nThe offset for the spline interpolation of the CEST data is:%4.2f Hz.\n',Minimum);

%-------------------------------------------------------------------------------------------
%Adjusted frequency range (shifted according to B0-correction)
Offset=-Minimum;
freqRange=linspace(-Range+Offset,Range+Offset,size(CESTdata_2,2)); 
 
freqRangeHR=linspace(-Range+Offset,Range+Offset,(2*Range+1)); 

pros = spline(freqRange,CESTdata_2,freqRangeHR); 

[min_val,idx]=min(pros(:)); 
[row,col]=ind2sub(size(pros),idx); 


%--------------------------------------------------------------

if(Seperate==1)
% Separate into 2 subspectra to correct for oscillations (=DOSE filtering)
ROI_Zspectrum = CESTdata_2;

l=1;
for i=1:2:length(frequenzen)
    freq_odd(1,l)=frequenzen(1,i)-Minimum;
    Z_odd(1,l)=ROI_Zspectrum(1,i);
    l=l+1;
end

%Add global minimum to both subspectra
stepsize=2*Range/(length(ROI_Zspectrum)-1);
if(abs(Minimum)~=stepsize && Minimum~=0 && abs(Minimum)~=2*stepsize && abs(Minimum)~=3*stepsize && abs(Minimum)~=4*stepsize && abs(Minimum)~=5*stepsize && abs(Minimum)~=6*stepsize && abs(Minimum)~=8*stepsize && abs(Minimum)~=12*stepsize && abs(Minimum)~=Range) %sonst doppelt 0 bei Offsets
freq_odd(1,length(freq_odd)+1)=0;
freq_odd=sort(freq_odd);
[mV,closestIndex_odd] = min(abs(0-freq_odd(1,:)));
Z_odd(1,(closestIndex_odd+1):length(Z_odd)+1)=Z_odd(1,closestIndex_odd:end); 
Z_odd(1,closestIndex_odd)=pros(1,col); 
end

l=1;
for i=2:2:length(frequenzen)-1
    freq_even(1,l)=frequenzen(1,i)-Minimum;
    Z_even(1,l)=ROI_Zspectrum(1,i);
    l=l+1;
end

%Add global minimum to both subspectra
if(Minimum~=0 && abs(Minimum)~=stepsize && abs(Minimum)~=2*stepsize && abs(Minimum)~=3*stepsize && abs(Minimum)~=4*stepsize && abs(Minimum)~=5*stepsize && abs(Minimum)~=7*stepsize && abs(Minimum)~=8*stepsize && abs(Minimum)~=9*stepsize && abs(Minimum)~=12*stepsize && abs(Minimum)~=Range) %sonst doppelt 0 bei Offsets
freq_even(1,length(freq_even)+1)=0;
freq_even=sort(freq_even);
[mV2,closestIndex_even] = min(abs(0-freq_even(1,:)));
Z_even(1,(closestIndex_even+1):length(Z_even)+1)=Z_even(1,closestIndex_even:end); 
Z_even(1,closestIndex_even)=pros(1,col); 
end


freq_spline_odd=linspace(min(freq_odd), max(freq_odd), max(freq_odd)-min(freq_odd)+1);
freq_spline_even=linspace(min(freq_even), max(freq_even), max(freq_even)-min(freq_even)+1);

pros_odd = spline(freq_odd,Z_odd,freq_spline_odd);
pros_even = spline(freq_even,Z_even,freq_spline_even);


%Determine the common frequency range for both subspectra
min_Frequenz(1,1)=freq_spline_odd(1,1); 
max_Frequenz(1,1)=freq_spline_odd(1,end); 
min_Frequenz(2,1)=freq_spline_even(1,1); 
max_Frequenz(2,1)=freq_spline_even(1,end); 

Untere_Frequenzgrenze=max(min_Frequenz);
Obere_Frequenzgrenze=min(max_Frequenz);


%Determine corresponding signal values
f_odd=freq_spline_odd(1,find(freq_spline_odd(1,:) == Untere_Frequenzgrenze):find(freq_spline_odd(1,:) == Obere_Frequenzgrenze)); 
Z_odd=pros_odd(1,find(freq_spline_odd(1,:) == Untere_Frequenzgrenze):find(freq_spline_odd(1,:) == Obere_Frequenzgrenze)); 

f_even=freq_spline_even(1,find(freq_spline_even(1,:) == Untere_Frequenzgrenze):find(freq_spline_even(1,:) == Obere_Frequenzgrenze)); 
Z_even=pros_even(1,find(freq_spline_even(1,:) == Untere_Frequenzgrenze):find(freq_spline_even(1,:) == Obere_Frequenzgrenze)); 


%Average both subspectra
for i=1:length(Z_even)
   mean_pros(1,i)=(Z_odd(1,i)+Z_even(1,i))/2; 
    
end


%B0_correction for the newly calculated averaged CEST spectrum
[min_val_mean,idx_mean]=min(mean_pros(:));
[row_mean,col_mean]=ind2sub(size(mean_pros),idx_mean);
Minimum_mean=(f_even(1,1)-1)+col_mean;
f_even=f_even-Minimum_mean;

Center=abs(Untere_Frequenzgrenze)+1+Minimum_mean;

%Save (pixelwise x-values = frequency infos; y-values = signal values; B0-offset)
Conclusion_seperate_x(j,k,1:length(f_even))=squeeze(f_even);
Conclusion_seperate_y(j,k,1:length(mean_pros))=squeeze(mean_pros);
Conclusion_Center(j,k,:)=Center;
end
%--------------------------------------------------------------

%Calculate pixelweise asymmetry spectra

None=length(freqRangeHR)/2; 

if (col > None && col ~= 2*None) 

for i=1:((2*None-1)-(col-1)) 
   asym(i)= pros(col-i)-pros(col+i);
end

freqRangeVER=linspace(0,Range-abs(Offset),(2*None-1)-(col-1)); 
[m, ind] = min(abs(freqRangeVER(1,:) - Int_Offset)); 
max_val=asym(1, ind);
end

if (col <= None && col ~=1) 
    
for i=1:((col-1)) 
   asym(i)= pros(col-i)-pros(col+i);
end

freqRangeVER=linspace(0,Range-abs(Offset),(col-1)); 
[m, ind] = min(abs(freqRangeVER(1,:) - Int_Offset)); 
max_val=asym(1, ind);
end

if (col==1) %Handle outliers (if the global minimum of a CEST spectrum was found at the outer sides of the frequency range)
AUC_12_22=0;
max_val=0;
asym=zeros(1,2000);
freqRangeVER=zeros(1,2000);
end

if (col==2*None) %Handle outliers (if the global minimum of a CEST spectrum was found at the outer sides of the frequency range)
AUC_12_22=0;
max_val=0;
asym=zeros(1,2000);
freqRangeVER=zeros(1,2000);
end

Help_a=1; 

%Only count pixels that completely contain the frequency range to be quantified (depending on pixelwise B0 correction)
if (length(freqRangeVER)>highBound_AUC && col ~= 2*None && col ~= 1)
[minValue_AUC3,closestIndex_AUC3] = min(abs(lowBound_AUC-freqRangeVER(1,:)));
[minValue_AUC4,closestIndex_AUC4] = min(abs(highBound_AUC-freqRangeVER(1,:)));
    if (closestIndex_AUC3 ~= closestIndex_AUC4) 
    AUC_12_22=(trapz(freqRangeVER(1,closestIndex_AUC3:closestIndex_AUC4),asym(1,closestIndex_AUC3:closestIndex_AUC4))); 
    Help_a=2;
    end
end

if (Help_a==1) 
AUC_12_22=0;
end

if (length(freqRangeVER)<=highBound_AUC && col ~= 1)    
AUC_12_22=0;
end

%---------------------------------------------

%Save all results (pixelwise)
Conclusion_frequencies(j,k,:)=squeeze(freqRange); %save CEST spectra (frequency infos) pixelwise
Conclusion_data(j,k,:)=squeeze(CESTdata_2); %save CEST spectra (signal values) pixelwiese

Conclusion_asym(j,k,:)=max_val; %save asym values pixelwise
Conclusion_AUC(j,k,:)=AUC_12_22; %save AUC values pixelwise
Conclusion_asymKurven(j,k,1:length(asym))=squeeze(asym); %save asym spectra values pixelweise
Conclusion_asymfrequencies(j,k,1:length(asym))=squeeze(freqRangeVER); %save asym spectra frequencies pixelwise

Conclusion_Minimum(j,k,:)=Minimum; %save B0-offset pixelwise
Conclusion_spline_x(j,k,:)=squeeze(freqRangeHR); %save splined CEST spectra (frequency infos) pixelwise
Conclusion_spline_y(j,k,:)=squeeze(pros); %save splined CEST spectra (signal values) pixelwise

Conclusion_sat(j,k,:)=min_val; %save global minimum signal value (saturation efficency)

end
end
end


%%

% Add 1: Calculation of global CEST spectrum for a selected layer and ROI (unfiltered)
%
% Be careful: As the pixelwise averaging is based on finding the common frequency range, pixels with significantly higher B0-offsets than the rest should be
% avoided (you can see the B0-offset distribution at "Conclusion_Minimum" and exclude individual pixels from the ROI (ROI(x,y)=0) if necessary).

for i=1:Pixelx
    for j=1:Pixely
    if (Conclusion_Minimum(i,j)==0)
       Conclusion_Minimum(i,j)=0.001; %This means that after the next calculation, pixels that are in the ROI can still be distinguished because they have the value 0.001 instead of 0!
    end 
    Conclusion_Minimum2(i,j)=Conclusion_Minimum(i,j).*ROI(i,j);
    end
end

%Determine the common frequency range for all pixels within the ROI
Conclusion_Minimum2(Conclusion_Minimum2==0)=NaN; 
LeftBound=min(min(Conclusion_Minimum2(:,:)));
RightBound=max(max(Conclusion_Minimum2(:,:)));

if(LeftBound>0) 
    LeftBound=0;
end
if(RightBound==0.001) 
    RightBound=0;
end

Span=2*Range-(abs(LeftBound)+abs(RightBound))+1;

freqRangesum2=linspace((-Range+abs(LeftBound)),(Range-abs(RightBound)),Span); %common frequency range


%Determine the corresponding signal values 
Conclusion_ver2=zeros(Pixelx,Pixely,Span);
clear i j k
for i=1:Pixelx
    for j=1:Pixely
        if(ROI(i,j) == 1)
        h1=squeeze(Conclusion_spline_y(i,j,:));
        h2=squeeze(Conclusion_spline_x(i,j,:));
        Conclusion_ver2(i,j,:)=squeeze(h1(find(h2 == (-Range+abs(LeftBound))):find(h2 == (Range-abs(RightBound)))));
        
        clear h1 h2
        end
    end
end


%Average signal values for each frequency offset
number = sum(ROI(:)~=0);
Z_averaged2=zeros(1,Span);

clear i j k
for k=1:Span
for i=1:Pixelx
    for j=1:Pixely
        
        if(ROI(i,j) == 1)
        Z_averaged2(1,k) = Conclusion_ver2(i,j,k)+Z_averaged2(1,k);     
        end  
        
    end
end
end

%Normalization based on ROI size
Z_averaged_norm2=Z_averaged2/number;

%Plot CEST spectrum (global; unfiltered)
figure();
plot(freqRangesum2/400,Z_averaged_norm2,'Linewidth',3)
xlabel('Frequency offset /ppm')
ylabel('Normalized signal')
set(gca,'FontSize',16)
xlim([-4 4])
ylim([0.0 1.0])
title('Global CEST spectrum (unfiltered)')



%%
% Add 2: Pixel-wise CEST contrast (based on unfiltered data)
% Idea: ROI serves as a mask here, therefore multiplication (component-wise) with ratio

ROI_copy=ROI;

%Sort out clusters from noise:
binaryImage = bwareaopen(Conclusion_asym, MinCluster);
ZZ=Conclusion_asym.*binaryImage;

for i=1:Pixelx
    for j=1:Pixely
        
       if (ROI(i,j)==1)
        ROI_data(i,j)=AA(i,j);  
        ROI_asym(i,j)=Conclusion_asym(i,j);
        ROI_AUC(i,j)=Conclusion_AUC(i,j);
       end
       
       if (ROI(i,j)==0)
       ROI_data(i,j,:)=AA(i,j)*0;  
       ROI_asym(i,j)=Conclusion_asym(i,j)*0.00000;
       ROI_AUC(i,j)=Conclusion_AUC(i,j)*0.00000;
       end

    end
end

%%Plot CEST maps (unfiltered; uncomment if necessary)
% figure;
% ax1 = axes;
% imagesc(ROI_data);
% colormap(ax1,'gray');
% 
% %single offset
% figure;
% ax1 = axes;
% imagesc(AA);
% title(sprintf('Asymmety value at %.0f Hz.', Int_Offset));
% colormap(ax1,'gray');
% ax2 = axes;
% imagesc(ax2,ZZ,'alphadata',ROI_asym~=0); 
% colormap(ax2,'jet');
% colormap(ax2,'jet');
% caxis(ax2,[min(nonzeros(ZZ)) max(nonzeros(ZZ))]);
% ax2.Visible = 'off';
% linkprop([ax1 ax2],'Position');
% colorbar;
% %caxis([min(ROI_asym(ROI_asym>0))-0.001 max(ROI_asym(ROI_asym>0))+0.001])
% %caxis([-0.05 0.15])
% 
% %AUC-based
% figure;
% ax1 = axes;
% imagesc(AA);
% title(sprintf('Area under curve at the frequency range %.0f Hz to %.0f Hz.', lowBound_AUC, highBound_AUC));
% colormap(ax1,'gray');
% ax2 = axes;
% imagesc(ax2,ROI_AUC,'alphadata',ROI_AUC~=0);
% colormap(ax2,'jet');
% colormap(ax2,'jet');
% caxis(ax2,[min(nonzeros(ROI_AUC)) max(nonzeros(ROI_AUC))]);
% ax2.Visible = 'off';
% linkprop([ax1 ax2],'Position');
% colorbar;
% %caxis([min(ROI_AUC(ROI_AUC>0))-1 max(ROI_AUC(ROI_AUC>0))+1])
% %caxis([-20 70])


%%

%Add 3: Determine global asymmetry spectrum (unfiltered) by averaging of all asymmetry curves within the selected ROI 

clear i j k

for i=1:Pixelx
    for j=1:Pixely
        if(ROI(i,j)==1)  
        Asym_data2(i,j,:)=Conclusion_asymKurven(i,j,:);
        Asym_Offset(i,j)=Conclusion_Minimum(i,j);
        end
        if(ROI(i,j)==0)
        Asym_data2(i,j,:)=Conclusion_asymKurven(i,j,:)*0;
        Asym_Offset(i,j)=0;
        end
    end

end


number = sum(ROI(:)~=0); %number of relevant pixels within ROI
ROI_Asymspectrum=zeros(1,Range);

AbsoluteOffsets=abs(Asym_Offset); %It only makes sense to calculate the mean value in the frequency range that all frequency lists of the asymmetry curves cover. For this reason, only frequencies from 0 to the maximum offset are included in ROI!
Bound=max(AbsoluteOffsets(:));

for k=1:(Range-0) %Calculate mean value within the ROI (sum over all elements, diveded by the number of elements which are not 0)
 for i=1:Pixelx
    for j=1:Pixely
    ROI_Asymspectrum(1,k)=Asym_data2(i,j,k)+squeeze(ROI_Asymspectrum(1,k));
    end
 end  
   ababa=(Asym_data2(:,:,k)); %number of contributions which contain a value unequal 0 for a certain frequency offset k
   No_unequal0(1,k)=nnz(ababa); 
   ROI_Asymspectrum(k)=ROI_Asymspectrum(k)/No_unequal0(1,k); 
    
end


%Plot asymmetry spectrum (global; unfiltered)
x_ROI_Asymspectrum=linspace(0,length(ROI_Asymspectrum),Range);
figure(); 
plot(x_ROI_Asymspectrum/400,ROI_Asymspectrum, 'Linewidth',3)
xlim([0 (Range-400)/400]);
xlabel('Frequency offset / ppm') 
ylabel('Asymmetry') 
set(gca,'FontSize',16)
title('Global Asymmetry spectrum (unfiltered)')

if (PraeMeas==0)
copy_ROI_Asymspectrum=ROI_Asymspectrum;
end

if (PraeMeas~=0)
figure(); 
plot(copy_ROI_Asymspectrum)
hold on
plot(ROI_Asymspectrum)
legend({'pre infusion', 'post infusion'},'Location','northeast') %,
xlim([0 (Range-400)]);
set(gca,'FontSize',16)
title ('Global asymmetry spectrum')
xlabel('Frequency offset / ppm') 
ylabel('Asymmetry') 
hold off
end


%---------------------------------
%Quantification of CEST contrast based on global asymmetry:

freq_asym=linspace(0,Range,Range+1); 
Bound =Int_Offset;
[minValue,closestIndex] = min(abs(Bound-freq_asym(1,:)));
Asym_512=ROI_Asymspectrum(1,closestIndex); %single offset value

%AuC range 1 (e.g. 0.5-2ppm):
Bound_AuC1=200;
Bound_AuC2=800;
[minValue1,closestIndex1] = min(abs(Bound_AuC1-freq_asym(1,:)));
[minValue2,closestIndex2] = min(abs(Bound_AuC2-freq_asym(1,:)));

%AuC range 2 (e.g. 1.2-2.2ppm):
Bound_AuC3=480;
Bound_AuC4=880;
[minValue3,closestIndex3] = min(abs(Bound_AuC3-freq_asym(1,:)));
[minValue4,closestIndex4] = min(abs(Bound_AuC4-freq_asym(1,:)));

%Calculate area under curve
Int_05_2=(trapz(freq_asym(1,closestIndex1:closestIndex2),ROI_Asymspectrum(1,closestIndex1:closestIndex2))); %0.5-2ppm 
Int_12_22=(trapz(freq_asym(1,closestIndex3:closestIndex4),ROI_Asymspectrum(1,closestIndex3:closestIndex4))); %1.2-2.2ppm 

%Conclusion
A_Ergebnisse_asym(1,1)=Asym_512;
A_Ergebnisse_asym(1,2)=Int_05_2/(Bound_AuC2-Bound_AuC1);
A_Ergebnisse_asym(1,3)=Int_12_22/(Bound_AuC4-Bound_AuC3);


%Save results
if (PraeMeas==0)
A_Prae_Ergebnisse=A_Ergebnisse_asym;
end

if (PraeMeas~=0)
A_Post_Ergebnisse=A_Ergebnisse_asym;
end


%%
%Add 4a: Processing of DOSE-filtered CEST spectra
%--> Pre-calculations & Pixelmaps

Conclusion_asym_mean=zeros(Pixelx,Pixely,Range);
asym_mean_512=zeros(Pixelx,Pixely,1);
Range_mean=Range-stepsize;
AUC_DOSE=zeros(Pixelx,Pixely);

%Calculation of asymmetry spectra pixelwise
for  k=1:Pixely
for  j=1:Pixelx
if(Conclusion_Center(j,k)<2*Range_mean && Conclusion_Center(j,k)>1) 
if(ROI(j,k)==1)
    Center=Conclusion_Center(j,k);  

if(Center<Range_mean)
Center2=Center; 
clear i asym_mean

for i=1:Center2-1 
   asym_mean(1,i)=squeeze(Conclusion_seperate_y(j,k,Center2-i))-squeeze(Conclusion_seperate_y(j,k,Center2+i));
end  
end


if(Center>=Range_mean) 
Center2=2*Range_mean-Center;
clear i asym_mean
for i=1:Center2-1 
   asym_mean(1,i)=squeeze(Conclusion_seperate_y(j,k,Center-i))-squeeze(Conclusion_seperate_y(j,k,Center+i)); 
end  

end
  
Conclusion_asym_mean(j,k,1:length(asym_mean))=squeeze(asym_mean);
asym_mean_512(j,k)=squeeze(Conclusion_asym_mean(j,k,Int_Offset));

freq_asym_help=linspace(1,length(asym_mean),length(asym_mean));
AUC_12_22=(trapz(freq_asym_help(1,lowBound_AUC:highBound_AUC),asym_mean(1,lowBound_AUC:highBound_AUC)))/(highBound_AUC-lowBound_AUC); %1.2-2.2ppm
clear freq_asym_help

AUC_DOSE(j,k)=AUC_12_22;
end
end
end
end

%CEST maps
binaryImage2 = bwareaopen(asym_mean_512, MinCluster);
ZZ2=asym_mean_512.*binaryImage2;
binaryImage3 = bwareaopen(AUC_DOSE, MinCluster);
ZZ3=AUC_DOSE.*binaryImage3;

%single offset
figure;
ax1 = axes;
imagesc(AA);
colormap(ax1,'gray');
ax2 = axes;
imagesc(ax2,ZZ2,'alphadata',asym_mean_512~=0); 
colormap(ax2,'jet');
colormap(ax2,'jet');
caxis(ax2,[min(nonzeros(ZZ2)) max(nonzeros(ZZ2))]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
colorbar;
caxis([-0.05 0.15])


%AUC
figure;
ax1 = axes;
imagesc(AA);
colormap(ax1,'gray');
ax2 = axes;
imagesc(ax2,ZZ3,'alphadata',AUC_DOSE~=0); 
colormap(ax2,'jet');
colormap(ax2,'jet');
caxis(ax2,[min(nonzeros(ZZ3)) max(nonzeros(ZZ3))]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
colorbar;
caxis([-0.1 0.15])

%Save
saveas(gcf,'Prae_AUC_Map.fig')
save("Prae_AUC_Map_Data","AUC_DOSE")


%--------------------
% %Averaging of pixelwise asymmetry spectra (similar to Add 3) 
% 
% clear i j k
% 
% for i=1:Pixelx
%     for j=1:Pixely
%         if(ROI(i,j)==1)  
%         Asym_mean_data(i,j,:)=Conclusion_asym_mean(i,j,:);
%         end
%         if(ROI(i,j)==0)
%         Asym_mean_data(i,j,:)=Conclusion_asymKurven(i,j,:)*0;
%         end
%     end
% 
% end
% 
% ROI_Asymspectrum_mean=zeros(1,Range_mean);
% 
% for k=1:Range_mean %Calculate mean value within the ROI (sum over all elements, diveded by the number of elements which are not 0)
%  for i=1:Pixelx
%     for j=1:Pixely
%     ROI_Asymspectrum_mean(1,k)=Asym_mean_data(i,j,k)+squeeze(ROI_Asymspectrum_mean(1,k));
%     end
%  end  
%    ababa_mean=(Asym_mean_data(:,:,k)); %number of contributions which contain a value unequal 0 for a certain frequency offset k
%    No_unequal0_mean(1,k)=nnz(ababa_mean); 
%    ROI_Asymspectrum_mean(k)=ROI_Asymspectrum_mean(k)/No_unequal0_mean(1,k);   
%     
% end
%
% figure(); 
% plot(ROI_Asymspectrum_mean)
% xlim([0 (Range_mean-400)]);
% %ylim([-0.03 0.08])
% title ('Asymmetriespektrum aus einzelnen Pixeln im gesamten ROI bei Seperate')
% xlabel('Frequenzoffset (Hz)') 
% ylabel('MTR_{Asym}') 


%%
%Add 4b: Processing the DOSE-filtered CEST spectra
%--> global CEST spectrum (DOSE-filtered)
%--> global asymmetry spectrum (DOSE-filtered)
%--> quantification


%Calculation of global CEST spectrum for a selected layer and ROI (filtered)
clear i j k

cop_Min=Conclusion_Minimum;
Conclusion_Minimum=Conclusion_Center-(Range_mean+1);

%Determine the common frequency range for all pixels within the ROI
for i=1:Pixelx
    for j=1:Pixely
    if (Conclusion_Minimum(i,j)==0)
       Conclusion_Minimum(i,j)=0.001; %This means that after the next calculation, pixels that are in the ROI can still be distinguished because they have the value 0.001 instead of 0!
    end 
    Conclusion_Minimum2(i,j)=Conclusion_Minimum(i,j).*ROI(i,j);
    end
end

Conclusion_Minimum2(Conclusion_Minimum2==0)=NaN; 
LeftBound=min(min(Conclusion_Minimum2(:,:)));
RightBound=max(max(Conclusion_Minimum2(:,:)));
Span=2*Range_mean-(abs(LeftBound)+abs(RightBound))+1;
freqRangesum=linspace((-Range_mean+abs(LeftBound)),(Range_mean-abs(RightBound)),Span);


%Determine the corresponding signal values 
Conclusion_ver=zeros(Pixelx,Pixely,Span);
clear i j k
for i=1:Pixelx
    for j=1:Pixely
        if(ROI(i,j) == 1)
        h1=squeeze(Conclusion_seperate_y(i,j,:));
        h2=squeeze(Conclusion_seperate_x(i,j,:));
        Conclusion_ver(i,j,:)=squeeze(h1(find(h2 == (-Range_mean+abs(LeftBound))):find(h2 == (Range_mean-abs(RightBound)))));
        clear h1 h2
        end
    end
end

%Average signal values for each frequency offset
number = sum(ROI(:)~=0);
Z_averaged=zeros(1,Span);

clear i j k
for k=1:Span
for i=1:Pixelx
    for j=1:Pixely
        
        if(ROI(i,j) == 1)
        Z_averaged(1,k) = Conclusion_ver(i,j,k)+Z_averaged(1,k);     
        end  
        
    end
end
end

%Save
save("Prae_y.mat","Conclusion_ver")
save("Prae_x.mat","freqRangesum")
save("Prae_ROI_complete.mat","ROI")

%Normalization based on ROI size
Z_averaged_norm=Z_averaged/number;

% %Plot CEST spectrum (global; filtered)
% figure();
% plot(freqRangesum,Z_averaged_norm,'Linewidth',3)
% saveas(gcf,'Prae_CEST-spectrum.fig')

%-------------------------------------------------------------------------
%Determine asymmetry
[min_val_new,idx_new]=min(Z_averaged_norm(:));
[row_new,col_new]=ind2sub(size(Z_averaged_norm),idx_new);

clear asym_new
for i=1:1000 %adjust asym calculation range if necessary
   asym_new(i)= Z_averaged_norm(col_new-i)-Z_averaged_norm(col_new+i);
end

% %Plot asymmetry spectrum (global; filtered)
% figure
% plot(asym_new)
% title('Asymmtrie aus DOSE (pixelweise)')
% xlabel('Frequenzoffset (Hz)') 
% ylabel('MTR_{Asym}') 
% saveas(gcf,'Prae_Asymmetriespektrum.fig')


%-------------------------------------------------------------------------
%Quantification
freq_asym=linspace(0,1100,1101); 

%Single offset
Bound =Int_Offset;
[minValue,closestIndex] = min(abs(Bound-freq_asym(1,:)));
Asym_512=asym_new(1,closestIndex); 

%AuC range 1 (e.g. 0.5-2ppm):
Bound_AuC1=200;
Bound_AuC2=800;
[minValue1,closestIndex1] = min(abs(Bound_AuC1-freq_asym(1,:)));
[minValue2,closestIndex2] = min(abs(Bound_AuC2-freq_asym(1,:)));

%AuC range 2 (e.g. 1.2-2.2ppm):
Bound_AuC3=480;
Bound_AuC4=880;
[minValue3,closestIndex3] = min(abs(Bound_AuC3-freq_asym(1,:)));
[minValue4,closestIndex4] = min(abs(Bound_AuC4-freq_asym(1,:)));

%Calculate AUC
Int_05_2=(trapz(freq_asym(1,closestIndex1:closestIndex2),asym_new(1,closestIndex1:closestIndex2))); 
Int_12_22=(trapz(freq_asym(1,closestIndex3:closestIndex4),asym_new(1,closestIndex3:closestIndex4))); 

%conclusion
A_DOSE_Ergebnisse_asym(1,1)=Asym_512;
A_DOSE_Ergebnisse_asym(2,1)=Int_05_2/(Bound_AuC2-Bound_AuC1);
A_DOSE_Ergebnisse_asym(3,1)=Int_12_22/(Bound_AuC4-Bound_AuC3);


%%
%Add 5a: Determination of outliers after DOSE filtering
%Based on marking pixels whose signal values are too extreme (negative or significantly above 1; values can be adjusted within the for-loops)

number_p=number;
OutLier=0;
OutLier2=0;
OutLier3=0;
OutLier_Position=zeros(Pixelx,Pixely);

for i=1:96
    for j=1:96
        if ROI(i,j)==1
            controlle_extrem=max(squeeze(Conclusion_ver(i,j,:)))
            if (max(squeeze(Conclusion_ver(i,j,:))>1.05) ||min(squeeze(Conclusion_ver(i,j,:))<-0.02))
               OutLier=OutLier+1; 
               OutLier_Position(i,j)=1;
            end
            if (max(squeeze(Conclusion_ver(i,j,:))>1.15) ||min(squeeze(Conclusion_ver(i,j,:))<-0.05))
               OutLier2=OutLier2+1; 
               OutLier_Position(i,j)=2;
            end
        end
    end
end

%Plot map higlighting outlier
figure;
ax1 = axes;
imagesc(AA);
colormap(ax1,'gray');
ax2 = axes;
imagesc(ax2,OutLier_Position,'alphadata',OutLier_Position~=0);
colormap(ax2,'jet');
colormap(ax2,'jet');
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
%colorbar;
caxis([0 2])
saveas(gcf,'Prae_Ausreiszer.fig')

Ausreisser_Stufe1=OutLier/number_p;
Ausreisser_Stufe2=OutLier2/number_p;


%Generate new ROI without these outliers
ROI_clean=ROI;
for i=1:Pixelx
    for j=1:Pixely
        if(OutLier_Position(i,j)==2)
            ROI_clean(i,j)=0;
        end
    end
end


%%
% Add 5b: Like 4b, but outliers are removed for calculations (= use new ROI "ROI_clean")
%--> global CEST spectrum (DOSE-filtered)
%--> global asymmetry spectrum (DOSE-filtered)
%--> quantification

%Calculation of global CEST spectrum for a selected layer and ROI (filtered; outliers excluded)
clear i j k
cop_Min=Conclusion_Minimum;
Conclusion_Minimum=Conclusion_Center-(Range_mean+1);

%Determine the common frequency range for all pixels within the ROI
for i=1:Pixelx
    for j=1:Pixely
    if (Conclusion_Minimum(i,j)==0)
       Conclusion_Minimum(i,j)=0.001; %This means that after the next calculation, pixels that are in the ROI can still be distinguished because they have the value 0.001 instead of 0!
    end 
    Conclusion_Minimum2(i,j)=Conclusion_Minimum(i,j).*ROI(i,j);
    end
end

Conclusion_Minimum2(Conclusion_Minimum2==0)=NaN; 
LeftBound=min(min(Conclusion_Minimum2(:,:)));
RightBound=max(max(Conclusion_Minimum2(:,:)));
Span=2*Range_mean-(abs(LeftBound)+abs(RightBound))+1;
freqRangesum=linspace((-Range_mean+abs(LeftBound)),(Range_mean-abs(RightBound)),Span); 


%Determine the corresponding signal values 
Conclusion_ver=zeros(Pixelx,Pixely,Span);
clear i j k
for i=1:Pixelx
    for j=1:Pixely
        if(ROI(i,j) == 1)
            if(OutLier_Position(i,j)~=2)
        h1=squeeze(Conclusion_seperate_y(i,j,:));
        h2=squeeze(Conclusion_seperate_x(i,j,:));
        Conclusion_ver(i,j,:)=squeeze(h1(find(h2 == (-Range_mean+abs(LeftBound))):find(h2 == (Range_mean-abs(RightBound)))));
        clear h1 h2
            end
        end
    end
end


%Average signal values for each frequency offset
number = sum(ROI(:)~=0);
Z_averaged=zeros(1,Span);

clear i j k
for k=1:Span
for i=1:Pixelx
    for j=1:Pixely
        
        if(ROI(i,j) == 1)
        Z_averaged(1,k) = Conclusion_ver(i,j,k)+Z_averaged(1,k);     
        end  
        
    end
end
end

%Normalization based on ROI size
Z_averaged_norm=Z_averaged/(number-OutLier2); 

%Plot CEST spectrum (global; filtered; outliers excluded)
figure(); 
plot(freqRangesum/400,Z_averaged_norm, 'Linewidth', 3)
xlabel('Frequency offset /ppm')
ylabel('Normalized signal')
set(gca,'FontSize',16)
xlim([-4 4])
ylim([0.0 1.0])
title('Global CEST spectrum (filtered; outlier excluded)')

%Save
save("FuerMittelung_x.mat", "freqRangesum")
save("FuerMittelung_y.mat", "Z_averaged_norm")


%Calculate Asymmetrie
[min_val_new,idx_new]=min(Z_averaged_norm(:));
[row_new,col_new]=ind2sub(size(Z_averaged_norm),idx_new);

clear asym_new
for i=1:1000 %adjust asym calculation range if necessary
   asym_new(i)= Z_averaged_norm(col_new-i)-Z_averaged_norm(col_new+i);
end

%Plot asymmetry spectrum (global; filtered; outliers excluded)
freq_help_asym=linspace(1,1000,1000);
figure
plot(freq_help_asym/400,asym_new,'Linewidth',3)
xlabel('Frequency offset / ppm') 
ylabel('Asymmetry') 
set(gca,'FontSize',16)
xlim([0 3])
ylim([-0.03 0.10])
title('Global asymmetry spectrum (filtered; outlier excluded)')

%---------------------------
%Quantification
freq_asym=linspace(0,1100,1101); 

%single offset
Bound =Int_Offset;
[minValue,closestIndex] = min(abs(Bound-freq_asym(1,:)));
Asym_512=asym_new(1,closestIndex); 

%AuC range 1 (e.g. 0.5-2ppm):
Bound_AuC1=200;
Bound_AuC2=800;
[minValue1,closestIndex1] = min(abs(Bound_AuC1-freq_asym(1,:)));
[minValue2,closestIndex2] = min(abs(Bound_AuC2-freq_asym(1,:)));

%AUC range 2 (e.g. 1.2-2.2ppm):
Bound_AuC3=480;
Bound_AuC4=880;
[minValue3,closestIndex3] = min(abs(Bound_AuC3-freq_asym(1,:)));
[minValue4,closestIndex4] = min(abs(Bound_AuC4-freq_asym(1,:)));

%Calculate area under curve: 
Int_05_2=(trapz(freq_asym(1,closestIndex1:closestIndex2),asym_new(1,closestIndex1:closestIndex2))); 
Int_12_22=(trapz(freq_asym(1,closestIndex3:closestIndex4),asym_new(1,closestIndex3:closestIndex4))); 

%Conclusion
A_cop=A_DOSE_Ergebnisse_asym;
A_DOSE_Ergebnisse_asym(1,1)=Asym_512;
A_DOSE_Ergebnisse_asym(2,1)=Int_05_2/(Bound_AuC2-Bound_AuC1);
A_DOSE_Ergebnisse_asym(3,1)=Int_12_22/(Bound_AuC4-Bound_AuC3);

return

