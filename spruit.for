C-----------------------------------------------------------------------
      SUBROUTINE SPRUIT(ITOP,DM,CAMBIOXI,teff,flux)!
C-----------------------------------------------------------------------
C Rutina para el cálculo el crecimiento del núcleo usando la cota
C superior de ingestión de He calculada por Spruit (2015).

C Antes llamar a la rutina CANTIDADES.

      IMPLICIT NONE

c     DEFINIR VARIABLES EN LOS COMMONS

      COMMON /B/ ESTR(NMESH,NCOL),EST1(NMESH,NCOL),EST2(NMESH,6),
     @           DIF(NMESH,4),U(4*NMESH),V(4*NMESH),W(4*NMESH),
     @           TOPE(2000,9),PSITOT(NMESH),DIF_FLUX,DIF_TEFF

      COMMON /CANT/ RHO_C(NMESH,3),DEL_C(NMESH,3),AD_C(NMESH,3),
     @              CP_C(NMESH,3),OP_C(NMESH,3),DGN(NMESH),S_C(NMESH),
     @              EN_C(NMESH,3),ENEU_C(NMESH,3),IEOS_N(NMESH),
     @              GRAD_RAD(NMESH),GRAD_TRUE(NMESH),PMOLEC(NMESH),
     @              BBETA(NMESH),DIF_COE_C(NMESH),GRAD_MU(NMESH)

C     Buscamos dónde comienza el borde del núcleo.
C     Para ello, primero leemos de EST1(I,6) si cada capa es radiativa
C     o convectiva. Para ello, usamos la siguiente convención:

C     1: transporte radiativo
C     2: convección
c     3: semiconveccion
c     4: mezcla termoalina
c     5: overshooting

      IEDGE1(I)= NINT(EST1(I,6))

C     Partimos de suponer que el borde está en la capa más externa

      n_borde = ITOP - 1

C     Vamos recorriendo las capas hasta encontrar la primera que es
C     convectiva

      DO WHILE (IEDGE1(n_borde) .EQ. 2)

            n_borde = n_borde - 1

      END DO

      n_borde = n_borde + 1  ! Para volver al lado convectivo del borde

C     Ahora que encontramos dónde está el borde, calculamos las
C     siguiente cantidades:

C     1) Luminosidad en el borde convectivo, en unidades de 10^33 erg/s.

      l_borde = (EST1(n_borde,1)+ESTR(n_borde,1)) * REFELE

C     2) Temperatura en el borde convectivo, en 10^6 K.

      T = exp(EST1(n_borde,4) + DLOG( 1.D0 + ESTR(n_borde,4) ))

C     3) Nabla adiabático y radiativo en el borde convectivo.

      nabla_ad = AD_C(n_borde,1)

      nabla_rad = GRAD_RAD(n_borde)

C     4) Calor específico en el borde convectivo.

      c_P = CP_C(n_borde,1)

C     5) Derivadas logarítmicas de la densidad respecto de T (delta)
C     y respecto de mu (phi)

      delta = DEL_C(n_borde,1)

      phi = 1

C     Abundancia de He en la capa radiativa inmediatamente afuera del
C     borde convectivo

      Y_afuera = ACMP(n_borde - 1, 2)

C     Promedio de la abundancia de He en el núcleo, pesado por la masa
C     en cada capa.

      sum = 0.d0

      DO i = n_borde, ITOP-1

            sum = sum + ( exp(PSITOT(i+1)) - exp(PSITOT(i)) ) * 
     @ ACMP(i, 2)

      END DO

      Y_core = (1.d0 - exp(PSITOT(n_borde)))^(-1) * sum

C     6) 

      mu_Y = 

C     Ahora que calculamos todas estas cantidades, evaluamos la
C     expresión de Spruit.

      Mdot = 0.5 * (delta * l_borde)/(phi * mu_Y * c_P * T) *
     @ (1 - nabla_ad/nabla_rad) / (Y_afuera - Y_core)

C     La masa a ingerir por el nucleo en un tiempo igual a STEP es:

      masa_a_ingerir = Mdot * STEP  ! UNIDADES???

C     Ahora que sabemos la masa ingerida, buscamos entre qué y qué capa
C     está contenida esa masa. Comenzamos suponiendo que solo se
C     ingiere la capa inmediatamente sobre el borde convectivo

      m_ingerida = (exp(PSITOT(n_borde)) - exp(PSITOT(n_borde-1))) * DM

      DO WHILE m_ingerida .lt. masa_a_ingerir

            n_borde = n_borde - 1

            m_ingerida = m_ingerida + (exp(PSITOT(n_borde)) - 
     @ exp(PSITOT(n_borde-1))) * DM
      END DO

C     Cuando acabe ese bucle, obtendremos en n_borde hasta qué capa
C     debería extenderse el borde convectivo. Luego, guardamos esta
C     información actualizada en EST1(I,6):

C     - De la capa 1 hasta la n_borde-1 la mezcla será adiabática:

      DO I = 1 , n_borde - 1

            EST1(I,6) = 1

      END DO

C     - De la capa n_borde hasta la ITOP-1 tendremos el nucleo
C     convectivo:


      DO I = n_borde, ITOP-1

            EST1(I,6) = 2

      END DO

      return
      END