%% Avaliação do Teste de Winston-Lutz
%Detecção automática Winston-Lutz
%Calibrado para WL com imagens do iView
%Campo quadrado (2.0cm x 2.0cm)
%Esfera de 8mm de diâmetro (Phantom Ball Bearing Elekta)

clear all
close all
DiametroEsfera_mm=8;
TamanhoCampo_mm=20;

% operador=input('Entre com o nome do Operador: ','s' );
data=date;
[Filename, Pathname] = uigetfile('*','Selecionar imagem do Winston-Lutz','Multiselect','on');


%%
if iscell(Filename)==1; %AJUSTAR PARA IMAGEM UNICA
    %%
   
    
    for i=1:length(Filename)
        Filename1=Filename{i};
        m=[Pathname Filename1];
        im = dicomread(m);
        info=dicominfo(m);
        orient=Filename{i};
        orient=orient(1:end-4);
        
        im2=im2double(im);
        size_pixel_mm=info.ImagePlanePixelSpacing(1);
        mag=info.RTImageSID/info.RadiationMachineSAD;
        
        size_field_px=TamanhoCampo_mm*mag/info.ImagePlanePixelSpacing(1);
        raio_esfera_px=((DiametroEsfera_mm/2)*mag)/info.ImagePlanePixelSpacing(1);
        
        campo=zeros(size(im2));
        campo(im2<0.64)=1;
        
        CC = bwconncomp(campo);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idx] = max(numPixels);
        campob=zeros(size(campo));
        campob(CC.PixelIdxList{idx}) = 1;
        campo=campob;
        
        esfera=ones(size(im2));
        esfera(im2<.45)=0;
        CC = bwconncomp(esfera);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idx] = max(numPixels);
        esfera(CC.PixelIdxList{idx}) = 0;
        
        
        %Centróide
        measurementscampo = regionprops(campo, 'centroid');
        centroidcampo = [measurementscampo.Centroid];
        centroidcampo=round(centroidcampo);
        
        measurementsesfera = regionprops(esfera, 'Centroid');
        centroidesfera = [measurementsesfera.Centroid];
        centroidesfera=round(centroidesfera);
        
        distancia_x_mm=(size_pixel_mm/mag)*(centroidcampo(1)-centroidesfera(1));
        distancia_y_mm=-(size_pixel_mm/mag)*(centroidcampo(2)-centroidesfera(2));
        distancia_mm= (size_pixel_mm/mag)*((((centroidcampo(1)-centroidesfera(1))^2)+((centroidcampo(2)-centroidesfera(2))^2))^.5);
        
        %%Ajuste manual
        
        %         resp3=input('Confirma detecção?: 0-não 1-OK ' );
        %
        %         if resp3==0
        %             figure,  imshow(im2,[]);
        %
        %             xlabel('Defina as bordas do campo manualmente');
        %             title([info.SeriesDescription(5:end)])
        %
        %
        %             mask_campo = imrect(gca,[512 512 size_field_px size_field_px]);
        %             wait(mask_campo);
        %             campo_manual=createMask(mask_campo);
        %
        %             xlabel('Defina a posição da esfera manualmente');
        %             mask_esfera =  imellipse(gca,[512 512 raio_esfera_px*2 raio_esfera_px*2]);
        %             wait(mask_esfera);
        %             esfera_manual=createMask(mask_esfera);
        %
        %
        %             measurementscampomanual = regionprops(campo_manual, 'Centroid');
        %             centroidcampomanual = [measurementscampomanual.Centroid];
        %             centroidcampo=round(centroidcampomanual);
        %             measurementsesferamanual = regionprops(esfera_manual, 'Centroid');
        %             centroidesferamanual = [measurementsesferamanual.Centroid];
        %             centroidesfera=round(centroidesferamanual);
        %     end
        %%%%%%%%%%%%%%%%%%
        
        
        ordem_imagens{i}=[orient];
        centroid_esfera_full(i,:)=centroidesfera;
        centroid_campo_full(i,:)=centroidcampo;
        imagens(:,:,i)=im2;
        distancia_mm_full(i)=distancia_mm;
        distancia_x_mm_full(i)=distancia_x_mm;
        distancia_y_mm_full(i)=distancia_y_mm;
    end
end

%%Imagem com resultados absolutos
figure
for i=1:length(Filename)

    subplot(3,3,i)
    imshow(imagens(:,:,i),[],'InitialMagnification',200);
    zoom(6);
    title(ordem_imagens{i});
    xlabel(distancia_x_mm_full(i));
    ylabel(distancia_y_mm_full(i));
    hold on
    plot(centroid_campo_full(i,1),centroid_campo_full(i,2),'xr');
    plot(centroid_esfera_full(i,1),centroid_esfera_full(i,2),'xg');
    viscircles(centroid_esfera_full(i,:),raio_esfera_px,'EdgeColor','g');
    rectangle('Position',[centroid_campo_full(i,1)-(size_field_px/2) centroid_campo_full(i,2)-(size_field_px/2) size_field_px size_field_px],'EdgeColor','r');
end


%Imprimir Tabela (Orientação)
tabela={'Campo', 'Long(mm)','Lat(mm)','Vert(mm)'};
for i=1:length(Filename)
    tabela{i+1,1}=ordem_imagens{i};
end

for i=1:length(Filename)
    if strcmp(tabela{i+1,1},'G0T0C0')==1
        tabela{i+1,2}=num2str(distancia_y_mm_full(i));
        tabela{i+1,3}=num2str(distancia_x_mm_full(i));
    else if strcmp(tabela{i+1,1},'G0T0C90')==1
            tabela{i+1,2}=num2str(distancia_y_mm_full(i));
            tabela{i+1,3}=num2str(distancia_x_mm_full(i));
        else if strcmp(tabela{i+1,1},'G0T0C270')==1
                tabela{i+1,2}=num2str(distancia_y_mm_full(i));
                tabela{i+1,3}=num2str(distancia_x_mm_full(i));
            else if strcmp(tabela{i+1,1},'G0T90')==1
                    tabela{i+1,3}=num2str(distancia_y_mm_full(i));
                    tabela{i+1,2}=num2str(-distancia_x_mm_full(i));
                else    if strcmp(tabela{i+1,1},'G0T270')==1
                        tabela{i+1,3}=num2str(-distancia_y_mm_full(i));
                        tabela{i+1,2}=num2str(distancia_x_mm_full(i));
                        
                    else if strcmp(tabela{i+1,1},'G90T0')==1
                            tabela{i+1,2}=num2str(distancia_y_mm_full(i));
                            tabela{i+1,4}=num2str(-distancia_x_mm_full(i));
                            
                        else if strcmp(tabela{i+1,1},'G180T0')==1
                                tabela{i+1,2}=num2str(distancia_y_mm_full(i));
                                tabela{i+1,3}=num2str(-distancia_x_mm_full(i));
                                
                            else if strcmp(tabela{i+1,1},'G270T0')==1
                                    tabela{i+1,2}=num2str(distancia_y_mm_full(i));
                                    tabela{i+1,4}=num2str(distancia_x_mm_full(i));
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end

end
tabela


save(['TesteWL_Boldrini_',info.StudyDate],'centroid_campo_full','Filename','centroid_esfera_full','distancia_x_mm_full','distancia_y_mm_full','DiametroEsfera_mm','imagens','mag','ordem_imagens','raio_esfera_px', 'size_field_px','size_pixel_mm','tabela','TamanhoCampo_mm')