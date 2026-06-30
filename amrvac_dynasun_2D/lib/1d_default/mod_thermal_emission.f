! module mod_thermal_emission -- synthesize emission flux of some
! thermal lines
! EUV lines database: 
! 'He_II_304' 'Fe_IX_171' 'Fe_XXIV_193' 'Fe_XIV_211' 'Fe_XVI_335'
! 'Fe_XVIII_94' 'Fe_XXI_131'
! subroutines: 
! get_EUV: get local EUV emission intensity (for 1d, 2d and 3d)
! get_SXR: get local Soft X-ray emission intensity (for 1d, 2d and 3d)

module mod_thermal_emission
  use mod_global_parameters
  use mod_physics

  implicit none

  integer :: n_aia
  double precision :: t_aia(1:101)
  double precision :: f_94(1:101),f_131(1:101),f_171(1:101)
  double precision :: f_193(1:101),f_211(1:101),f_304(1:101)
  double precision :: f_335(1:101)
  integer :: n_iris
  double precision :: t_iris(1:41)
  double precision :: f_1354(1:41)
  integer :: n_eis
  double precision :: t_eis1(1:60),t_eis2(1:60)
  double precision :: f_263(1:60),f_264(1:60),f_192(1:60),f_255(1:60)


  double precision :: vec_xI1(1:3),vec_xI2(1:3),vec_LOS(1:3)

  data n_aia / 101 /

  data t_aia / 4. ,  4.05, 4.1,  4.15, 4.2,  4.25, 4.3,  4.35, 4.4,  4.45, 4.5,&
       4.55, 4.6,  4.65, 4.7,  4.75, 4.8,  4.85, 4.9,  4.95, 5. ,  5.05, 5.1,&
       5.15, 5.2,  5.25, 5.3,  5.35, 5.4,  5.45, 5.5,  5.55, 5.6,  5.65, 5.7,&
       5.75, 5.8,  5.85, 5.9,  5.95, 6. ,  6.05, 6.1,  6.15, 6.2,  6.25, 6.3,&
       6.35, 6.4,  6.45, 6.5,  6.55, 6.6,  6.65, 6.7,  6.75, 6.8,  6.85, 6.9,&
       6.95, 7. ,  7.05, 7.1,  7.15, 7.2,  7.25, 7.3,  7.35, 7.4,  7.45, 7.5,&
       7.55, 7.6,  7.65, 7.7,  7.75, 7.8,  7.85, 7.9,  7.95, 8. ,  8.05, 8.1,&
       8.15, 8.2,  8.25, 8.3,  8.35, 8.4,  8.45, 8.5,  8.55, 8.6,  8.65, 8.7,&
       8.75, 8.8,  8.85, 8.9,  8.95, 9. /

  data f_94 / 4.25022959d-37, 4.35880298d-36, 3.57054296d-35, 2.18175426d-34,&
      8.97592571d-34, 2.68512961d-33, 7.49559346d-33, 2.11603751d-32,&
      5.39752853d-32, 1.02935904d-31, 1.33822307d-31, 1.40884290d-31,&
      1.54933156d-31, 2.07543102d-31, 3.42026227d-31, 6.31171444d-31,&
      1.16559416d-30, 1.95360497d-30, 2.77818735d-30, 3.43552578d-30,&
      4.04061803d-30, 4.75470982d-30, 5.65553769d-30, 6.70595782d-30,&
      7.80680354d-30, 8.93247715d-30, 1.02618156d-29, 1.25979030d-29,&
      1.88526483d-29, 3.62448572d-29, 7.50553279d-29, 1.42337571d-28,&
      2.37912813d-28, 3.55232305d-28, 4.84985757d-28, 6.20662827d-28,&
      7.66193687d-28, 9.30403645d-28, 1.10519802d-27, 1.25786927d-27,&
      1.34362634d-27, 1.33185242d-27, 1.22302081d-27, 1.05677973d-27,&
      9.23064720d-28, 8.78570994d-28, 8.02397416d-28, 5.87681142d-28,&
      3.82272695d-28, 3.11492649d-28, 3.85736090d-28, 5.98893519d-28,&
      9.57553548d-28, 1.46650267d-27, 2.10365847d-27, 2.79406671d-27,&
      3.39420087d-27, 3.71077520d-27, 3.57296767d-27, 2.95114380d-27,&
      2.02913103d-27, 1.13361825d-27, 5.13405629d-28, 2.01305089d-28,&
      8.15781482d-29, 4.28366817d-29, 3.08701543d-29, 2.68693906d-29,&
      2.51764203d-29, 2.41773103d-29, 2.33996083d-29, 2.26997246d-29,&
      2.20316143d-29, 2.13810001d-29, 2.07424438d-29, 2.01149189d-29,&
      1.94980213d-29, 1.88917920d-29, 1.82963583d-29, 1.77116920d-29,&
      1.71374392d-29, 1.65740593d-29, 1.60214447d-29, 1.54803205d-29,&
      1.49510777d-29, 1.44346818d-29, 1.39322305d-29, 1.34441897d-29,&
      1.29713709d-29, 1.25132618d-29, 1.20686068d-29, 1.14226584d-29,&
      1.09866413d-29, 1.05635524d-29, 1.01532444d-29, 9.75577134d-30,&
      9.37102736d-30, 8.99873335d-30, 8.63860172d-30, 8.29051944d-30,&
      7.95414793d-30 /

  data f_131 / 3.18403601d-37,   3.22254703d-36,   2.61657920d-35,&
      1.59575286d-34,   6.65779556d-34,   2.07015132d-33, 6.05768615d-33,&
        1.76074833d-32,   4.52633001d-32, 8.57121883d-32,   1.09184271d-31,&
        1.10207963d-31, 1.11371658d-31,   1.29105226d-31,   1.80385897d-31,&
      3.27295431d-31,   8.92002136d-31,   3.15214579d-30, 9.73440787d-30,&
        2.22709702d-29,   4.01788984d-29, 6.27471832d-29,   8.91764995d-29,&
        1.18725647d-28, 1.52888040d-28,   2.05082946d-28,   3.47651873d-28,&
      8.80482184d-28,   2.66533063d-27,   7.05805149d-27, 1.46072515d-26,&
        2.45282476d-26,   3.55303726d-26, 4.59075911d-26,   5.36503515d-26,&
        5.68444094d-26, 5.47222296d-26,   4.81119761d-26,   3.85959059d-26,&
      2.80383406d-26,   1.83977650d-26,   1.11182849d-26, 6.50748885d-27,&
        3.96843481d-27,   2.61876319d-27, 1.85525324d-27,   1.39717024d-27,&
        1.11504283d-27, 9.38169611d-28,   8.24801234d-28,   7.43331919d-28,&
      6.74537063d-28,   6.14495760d-28,   5.70805277d-28, 5.61219786d-28,&
        6.31981777d-28,   9.19747307d-28, 1.76795732d-27,   3.77985446d-27,&
        7.43166191d-27, 1.19785603d-26,   1.48234676d-26,   1.36673114d-26,&
      9.61047146d-27,   5.61209353d-27,   3.04779780d-27, 1.69378976d-27,&
        1.02113491d-27,   6.82223774d-28, 5.02099099d-28,   3.99377760d-28,&
        3.36279037d-28, 2.94767378d-28,   2.65740865d-28,   2.44396277d-28,&
      2.28003967d-28,   2.14941419d-28,   2.04178995d-28, 1.95031045d-28,&
        1.87011994d-28,   1.79777869d-28, 1.73093957d-28,   1.66795789d-28,&
        1.60785455d-28, 1.55002399d-28,   1.49418229d-28,   1.44022426d-28,&
      1.38807103d-28,   1.33772767d-28,   1.28908404d-28, 1.24196208d-28,&
        1.17437501d-28,   1.12854330d-28, 1.08410498d-28,   1.04112003d-28,&
        9.99529904d-29, 9.59358806d-29,   9.20512291d-29,   8.83009123d-29,&
      8.46817043d-29,   8.11921928d-29 /

  data f_171 / 2.98015581d-42, 1.24696230d-40, 3.37614652d-39, 5.64103034d-38,&
      5.20550266d-37, 2.77785939d-36, 1.16283616d-35, 6.50007689d-35,&
      9.96177399d-34, 1.89586076d-32, 2.10982799d-31, 1.36946479d-30,&
      6.27396553d-30, 2.29955134d-29, 7.13430211d-29, 1.91024282d-28,&
      4.35358848d-28, 7.94807808d-28, 1.07431875d-27, 1.08399488d-27,&
      9.16212938d-28, 7.34715770d-28, 6.59246382d-28, 9.13541375d-28,&
      2.05939035d-27, 5.08206555d-27, 1.10148083d-26, 2.01884662d-26,&
      3.13578384d-26, 4.14367719d-26, 5.36067711d-26, 8.74170213d-26,&
      1.64161233d-25, 2.94587860d-25, 4.76298332d-25, 6.91765639d-25,&
      9.08825111d-25, 1.08496183d-24, 1.17440114d-24, 1.13943939d-24,&
      9.71696981d-25, 7.09593688d-25, 4.31376399d-25, 2.12708486d-25,&
      8.47429567d-26, 3.17608104d-26, 1.95898842d-26, 1.98064242d-26,&
      1.67706555d-26, 8.99126003d-27, 3.29773878d-27, 1.28896127d-27,&
      8.51169698d-28, 7.53520167d-28, 6.18268143d-28, 4.30034650d-28,&
      2.78152409d-28, 1.95437088d-28, 1.65896278d-28, 1.68740181d-28,&
      1.76054383d-28, 1.63978419d-28, 1.32880591d-28, 1.00833205d-28,&
      7.82252806d-29, 6.36181741d-29, 5.34633869d-29, 4.58013864d-29,&
      3.97833422d-29, 3.49414760d-29, 3.09790940d-29, 2.76786227d-29,&
      2.48806269d-29, 2.24823367d-29, 2.04016653d-29, 1.85977413d-29,&
      1.70367499d-29, 1.56966125d-29, 1.45570643d-29, 1.35964565d-29,&
      1.27879263d-29, 1.21016980d-29, 1.15132499d-29, 1.09959628d-29,&
      1.05307482d-29, 1.01040261d-29, 9.70657096d-30, 9.33214234d-30,&
      8.97689427d-30, 8.63761192d-30, 8.31149879d-30, 7.85162401d-30,&
      7.53828281d-30, 7.23559452d-30, 6.94341530d-30, 6.66137038d-30,&
      6.38929156d-30, 6.12669083d-30, 5.87346434d-30, 5.62943622d-30,&
      5.39435202d-30 /

  data f_193 / 6.40066486d-32, 4.92737300d-31, 2.95342934d-30, 1.28061594d-29,&
      3.47747667d-29, 5.88554792d-29, 7.72171179d-29, 9.75609282d-29,&
      1.34318963d-28, 1.96252638d-28, 2.70163878d-28, 3.63192965d-28,&
      5.28087341d-28, 8.37821446d-28, 1.39089159d-27, 2.31749718d-27,&
      3.77510689d-27, 5.85198594d-27, 8.26021568d-27, 1.04870405d-26,&
      1.25209374d-26, 1.47406787d-26, 1.77174067d-26, 2.24098537d-26,&
      3.05926105d-26, 4.50018853d-26, 6.84720216d-26, 1.00595861d-25,&
      1.30759222d-25, 1.36481773d-25, 1.15943558d-25, 1.01467304d-25,&
      1.04092532d-25, 1.15071251d-25, 1.27416033d-25, 1.38463476d-25,&
      1.47882726d-25, 1.57041238d-25, 1.69786224d-25, 1.94970397d-25,&
      2.50332918d-25, 3.58321431d-25, 5.18061550d-25, 6.60405549d-25,&
      6.64085365d-25, 4.83825816d-25, 2.40545020d-25, 8.59534098d-26,&
      2.90920638d-26, 1.33204845d-26, 9.03933926d-27, 7.78910836d-27,&
      7.29342321d-27, 7.40267022d-27, 8.05279981d-27, 8.13829291d-27,&
      6.92634262d-27, 5.12521880d-27, 3.59527615d-27, 2.69617560d-27,&
      2.84432713d-27, 5.06697306d-27, 1.01281903d-26, 1.63526978d-26,&
      2.06759342d-26, 2.19482312d-26, 2.10050611d-26, 1.89837248d-26,&
      1.66347131d-26, 1.43071097d-26, 1.21518419d-26, 1.02078343d-26,&
      8.46936184d-27, 6.93015742d-27, 5.56973237d-27, 4.38951754d-27,&
      3.38456457d-27, 2.55309556d-27, 1.88904224d-27, 1.38057546d-27,&
      1.00718330d-27, 7.43581116d-28, 5.63562931d-28, 4.43359435d-28,&
      3.63923535d-28, 3.11248143d-28, 2.75586846d-28, 2.50672237d-28,&
      2.32419348d-28, 2.18325682d-28, 2.06834486d-28, 1.93497044d-28,&
      1.84540751d-28, 1.76356504d-28, 1.68741425d-28, 1.61566157d-28,&
      1.54754523d-28, 1.48249410d-28, 1.42020176d-28, 1.36045230d-28,&
      1.30307965d-28 /

  data f_211 / 4.74439912d-42, 1.95251522d-40, 5.19700194d-39, 8.53120166d-38,&
      7.72745727d-37, 4.04158559d-36, 1.64853511d-35, 8.56295439d-35,&
      1.17529722d-33, 2.16867729d-32, 2.40472264d-31, 1.56418133d-30,&
      7.20032889d-30, 2.65838271d-29, 8.33196904d-29, 2.26128236d-28,&
      5.24295811d-28, 9.77791121d-28, 1.35913489d-27, 1.43957785d-27,&
      1.37591544d-27, 1.49029886d-27, 2.06183401d-27, 3.31440622d-27,&
      5.42497318d-27, 8.41100374d-27, 1.17941366d-26, 1.49269794d-26,&
      1.71506074d-26, 1.71266353d-26, 1.51434781d-26, 1.36766622d-26,&
      1.33483562d-26, 1.36834518d-26, 1.45829002d-26, 1.62575306d-26,&
      1.88773347d-26, 2.22026986d-26, 2.54930499d-26, 2.80758138d-26,&
      3.06176409d-26, 3.62799792d-26, 5.13226109d-26, 8.46260744d-26,&
      1.38486586d-25, 1.86192535d-25, 1.78007934d-25, 1.16548409d-25,&
      5.89293257d-26, 2.69952884d-26, 1.24891081d-26, 6.41273176d-27,&
      4.08282914d-27, 3.26463328d-27, 2.76230280d-27, 2.08986882d-27,&
      1.37658470d-27, 8.48489381d-28, 5.19304217d-28, 3.19312514d-28,&
      2.02968197d-28, 1.50171666d-28, 1.39164218d-28, 1.42448821d-28,&
      1.41714519d-28, 1.33341059d-28, 1.20759270d-28, 1.07259692d-28,&
      9.44895400d-29, 8.29030041d-29, 7.25440631d-29, 6.33479483d-29,&
      5.51563757d-29, 4.79002469d-29, 4.14990482d-29, 3.59384972d-29,&
      3.12010860d-29, 2.72624742d-29, 2.40734791d-29, 2.15543565d-29,&
      1.95921688d-29, 1.80682882d-29, 1.68695662d-29, 1.59020936d-29,&
      1.50940886d-29, 1.43956179d-29, 1.37731622d-29, 1.32049043d-29,&
      1.26771875d-29, 1.21803879d-29, 1.17074716d-29, 1.10507836d-29,&
      1.06022834d-29, 1.01703080d-29, 9.75436986d-30, 9.35349257d-30,&
      8.96744546d-30, 8.59527489d-30, 8.23678940d-30, 7.89144480d-30,&
      7.55891138d-30 /

  data f_304 / 3.62695850d-32, 2.79969087d-31, 1.68340584d-30, 7.32681440d-30,&
      1.99967770d-29, 3.41296785d-29, 4.55409104d-29, 5.94994635d-29,&
      8.59864963d-29, 1.39787633d-28, 3.17701965d-28, 1.14474920d-27,&
      4.44845958d-27, 1.54785841d-26, 4.70265345d-26, 1.24524365d-25,&
      2.81535352d-25, 5.10093666d-25, 6.83545307d-25, 6.82110329d-25,&
      5.66886188d-25, 4.36205513d-25, 3.29265688d-25, 2.49802368d-25,&
      1.92527113d-25, 1.51058572d-25, 1.20596047d-25, 9.76884267d-26,&
      7.89979266d-26, 6.18224289d-26, 4.67298332d-26, 3.57934505d-26,&
      2.84535785d-26, 2.32853022d-26, 1.95228514d-26, 1.67880071d-26,&
      1.47608785d-26, 1.32199691d-26, 1.20070960d-26, 1.09378177d-26,&
      1.00031730d-26, 9.62434001d-27, 1.05063954d-26, 1.27267143d-26,&
      1.45923057d-26, 1.36746707d-26, 1.03466970d-26, 6.97647829d-27,&
      4.63141039d-27, 3.19031994d-27, 2.33373613d-27, 1.81589079d-27,&
      1.48446917d-27, 1.26611478d-27, 1.12617468d-27, 1.03625148d-27,&
      9.61400595d-28, 8.79016231d-28, 7.82612130d-28, 6.73762960d-28,&
      5.59717956d-28, 4.53010243d-28, 3.65712196d-28, 3.00958686d-28,&
      2.54011502d-28, 2.18102277d-28, 1.88736437d-28, 1.63817539d-28,&
      1.42283147d-28, 1.23631916d-28, 1.07526003d-28, 9.36797928d-29,&
      8.18565660d-29, 7.18152734d-29, 6.32523238d-29, 5.59513985d-29,&
      4.96614048d-29, 4.42518826d-29, 3.95487628d-29, 3.54690294d-29,&
      3.18953930d-29, 2.87720933d-29, 2.60186750d-29, 2.36011522d-29,&
      2.14717806d-29, 1.95905217d-29, 1.79287981d-29, 1.64562262d-29,&
      1.51489425d-29, 1.39876064d-29, 1.29496850d-29, 1.18665438d-29,&
      1.10240474d-29, 1.02643099d-29, 9.57780996d-30, 8.95465151d-30,&
      8.38950190d-30, 7.87283711d-30, 7.40136507d-30, 6.96804279d-30,&
      6.56945323d-30 /

  data f_335 / 2.46882661d-32, 1.89476632d-31, 1.13216502d-30, 4.89532008d-30,&
      1.32745970d-29, 2.25390335d-29, 3.00511672d-29, 3.96035934d-29,&
      5.77977656d-29, 8.58600736d-29, 1.14083000d-28, 1.48644411d-28,&
      2.15788823d-28, 3.51628877d-28, 6.12200698d-28, 1.08184987d-27,&
      1.85590697d-27, 2.91679107d-27, 3.94405396d-27, 4.63610680d-27,&
      5.13824456d-27, 5.66602209d-27, 6.30009232d-27, 7.03422868d-27,&
      7.77973918d-27, 8.32371831d-27, 8.56724316d-27, 8.62601374d-27,&
      8.13308844d-27, 6.53188216d-27, 4.55197029d-27, 3.57590087d-27,&
      3.59571707d-27, 4.03502770d-27, 4.54366411d-27, 4.96914990d-27,&
      5.24601170d-27, 5.39979250d-27, 5.43023669d-27, 5.26235042d-27,&
      4.91585495d-27, 4.52628362d-27, 4.13385020d-27, 3.67538967d-27,&
      3.39939742d-27, 3.81284533d-27, 5.02332701d-27, 6.19438602d-27,&
      6.49613071d-27, 6.04010475d-27, 5.24664275d-27, 4.37225997d-27,&
      3.52957182d-27, 2.76212276d-27, 2.08473158d-27, 1.50850518d-27,&
      1.04602472d-27, 7.13091243d-28, 5.34289645d-28, 5.21079581d-28,&
      6.22246365d-28, 6.99555864d-28, 6.29665489d-28, 4.45077026d-28,&
      2.67046793d-28, 1.52774686d-28, 9.18061770d-29, 6.09116074d-29,&
      4.48562572d-29, 3.59463696d-29, 3.05820218d-29, 2.70766652d-29,&
      2.46144034d-29, 2.27758450d-29, 2.13331183d-29, 2.01537836d-29,&
      1.91566180d-29, 1.82893912d-29, 1.75167748d-29, 1.68136168d-29,&
      1.61615595d-29, 1.55481846d-29, 1.49643236d-29, 1.44046656d-29,&
      1.38657085d-29, 1.33459068d-29, 1.28447380d-29, 1.23615682d-29,&
      1.18963296d-29, 1.14478976d-29, 1.10146637d-29, 1.04039479d-29,&
      9.98611410d-30, 9.58205147d-30, 9.19202009d-30, 8.81551313d-30,&
      8.45252127d-30, 8.10224764d-30, 7.76469090d-30, 7.43954323d-30,&
      7.12653873d-30 /


  data n_iris / 41 /

  data t_iris / 4.        , 4.1       , 4.2       , 4.3       , 4.40000001,&
      4.50000001, 4.60000001, 4.70000001, 4.80000001, 4.90000001, 5.00000001,&
      5.10000002, 5.20000002, 5.30000002, 5.40000002, 5.50000002, 5.60000002,&
      5.70000003, 5.80000003, 5.90000003, 6.00000003, 6.10000003, 6.20000003,&
      6.30000003, 6.40000004, 6.50000004, 6.60000004, 6.70000004, 6.80000004,&
      6.90000004, 7.00000004, 7.10000005, 7.20000005, 7.30000005, 7.40000005,&
      7.50000005, 7.60000005, 7.70000006, 7.80000006, 7.90000006, 8.00000006 /

  data f_1354 / 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 1.09503647d-39,&
      5.47214550d-36, 2.42433983d-33, 2.75295034d-31, 1.21929718d-29,&
      2.48392125d-28, 2.33268145d-27, 8.68623633d-27, 1.00166284d-26,&
      3.63126633d-27, 7.45174807d-28, 1.38224064d-28, 2.69270994d-29,&
      5.53314977d-30, 1.15313092d-30, 2.34195788d-31, 4.48242942d-32,&
      7.94976380d-33 /


  data n_eis  / 60 /

  data t_eis1 / 1.99526231d+05, 2.23872114d+05, 2.51188643d+05, 2.81838293d+05,&
      3.16227766d+05, 3.54813389d+05, 3.98107171d+05, 4.46683592d+05,&
      5.01187234d+05, 5.62341325d+05, 6.30957344d+05, 7.07945784d+05,&
      7.94328235d+05, 8.91250938d+05, 1.00000000d+06, 1.12201845d+06,&
      1.25892541d+06, 1.41253754d+06, 1.58489319d+06, 1.77827941d+06,&
      1.99526231d+06, 2.23872114d+06, 2.51188643d+06, 2.81838293d+06,&
      3.16227766d+06, 3.54813389d+06, 3.98107171d+06, 4.46683592d+06,&
      5.01187234d+06, 5.62341325d+06, 6.30957344d+06, 7.07945784d+06,&
      7.94328235d+06, 8.91250938d+06, 1.00000000d+07, 1.12201845d+07,&
      1.25892541d+07, 1.41253754d+07, 1.58489319d+07, 1.77827941d+07,&
      1.99526231d+07, 2.23872114d+07, 2.51188643d+07, 2.81838293d+07,&
      3.16227766d+07, 3.54813389d+07, 3.98107171d+07, 4.46683592d+07,&
      5.01187234d+07, 5.62341325d+07, 6.30957344d+07, 7.07945784d+07,&
      7.94328235d+07, 8.91250938d+07, 1.00000000d+08, 1.12201845d+08,&
      1.25892541d+08, 1.41253754d+08, 1.58489319d+08, 1.77827941d+08 /

  data t_eis2 / 1.99526231d+06, 2.23872114d+06, 2.51188643d+06, 2.81838293d+06,&
      3.16227766d+06, 3.54813389d+06, 3.98107171d+06, 4.46683592d+06,&
      5.01187234d+06, 5.62341325d+06, 6.30957344d+06, 7.07945784d+06,&
      7.94328235d+06, 8.91250938d+06, 1.00000000d+07, 1.12201845d+07,&
      1.25892541d+07, 1.41253754d+07, 1.58489319d+07, 1.77827941d+07,&
      1.99526231d+07, 2.23872114d+07, 2.51188643d+07, 2.81838293d+07,&
      3.16227766d+07, 3.54813389d+07, 3.98107171d+07, 4.46683592d+07,&
      5.01187234d+07, 5.62341325d+07, 6.30957344d+07, 7.07945784d+07,&
      7.94328235d+07, 8.91250938d+07, 1.00000000d+08, 1.12201845d+08,&
      1.25892541d+08, 1.41253754d+08, 1.58489319d+08, 1.77827941d+08,&
      1.99526231d+08, 2.23872114d+08, 2.51188643d+08, 2.81838293d+08,&
      3.16227766d+08, 3.54813389d+08, 3.98107171d+08, 4.46683592d+08,&
      5.01187234d+08, 5.62341325d+08, 6.30957344d+08, 7.07945784d+08,&
      7.94328235d+08, 8.91250938d+08, 1.00000000d+09, 1.12201845d+09,&
      1.25892541d+09, 1.41253754d+09, 1.58489319d+09, 1.77827941d+09 /

  data f_263 /  0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,&
      0.00000000d+00, 4.46454917d-45, 3.26774829d-42, 1.25292566d-39,&
      2.66922338d-37, 3.28497742d-35, 2.38677554d-33, 1.03937729d-31,&
      2.75075687d-30, 4.47961733d-29, 4.46729177d-28, 2.64862689d-27,&
      8.90863800d-27, 1.72437548d-26, 2.22217752d-26, 2.27999477d-26,&
      2.08264363d-26, 1.78226687d-26, 1.45821699d-26, 1.14675379d-26,&
      8.63082492d-27, 6.15925429d-27, 4.11252514d-27, 2.51530564d-27,&
      1.37090986d-27, 6.42443134d-28, 2.48392636d-28, 7.59187874d-29,&
      1.77852938d-29, 3.23945221d-30, 4.90533903d-31, 6.75458158d-32,&
      9.06878868d-33, 1.23927474d-33, 1.75769395d-34, 2.60710914d-35,&
      4.04318030d-36, 6.53500581d-37, 1.09365022d-37, 1.88383322d-38,&
      3.31425233d-39, 5.90964084d-40, 1.06147549d-40, 1.90706170d-41,&
      3.41331584d-42, 6.07310718d-43, 1.07364738d-43, 1.89085498d-44,&
      3.32598922d-45, 5.87125640d-46, 0.00000000d+00, 0.00000000d+00 /

  data f_264 /  0.00000000d+00, 2.81670057d-46, 1.28007268d-43, 2.54586603d-41,&
      2.67887256d-39, 1.68413285d-37, 6.85702304d-36, 1.91797284d-34,&
      3.84675839d-33, 5.69939170d-32, 6.36224608d-31, 5.39176489d-30,&
      3.45478458d-29, 1.64848693d-28, 5.71476364d-28, 1.39909997d-27,&
      2.37743056d-27, 2.86712530d-27, 2.65206348d-27, 2.07175767d-27,&
      1.47866767d-27, 1.01087374d-27, 6.79605811d-28, 4.54746770d-28,&
      3.04351751d-28, 2.03639149d-28, 1.35940991d-28, 9.01451939d-29,&
      5.91289972d-29, 3.81821178d-29, 2.41434696d-29, 1.48871220d-29,&
      8.93362094d-30, 5.21097445d-30, 2.95964719d-30, 1.64278748d-30,&
      8.95571660d-31, 4.82096011d-31, 2.57390991d-31, 1.36821781d-31,&
      7.27136350d-32, 3.87019426d-32, 2.06883430d-32, 1.11228884d-32,&
      6.01883313d-33, 3.27790676d-33, 1.79805012d-33, 9.93085346d-34,&
      5.52139556d-34, 3.08881387d-34, 1.73890315d-34, 9.84434964d-35,&
      5.60603378d-35, 3.20626492d-35, 1.84111068d-35, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00 /

  data f_192 /  0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 4.35772105d-44,&
      1.26162319d-41, 1.97471205d-39, 1.83409019d-37, 1.08206288d-35,&
      4.27914363d-34, 1.17943846d-32, 2.32565755d-31, 3.33087991d-30,&
      3.47013260d-29, 2.60375866d-28, 1.37737127d-27, 5.01053913d-27,&
      1.23479810d-26, 2.11310542d-26, 2.71831513d-26, 2.89851163d-26,&
      2.77312376d-26, 2.50025229d-26, 2.18323661d-26, 1.86980322d-26,&
      1.58035034d-26, 1.31985651d-26, 1.08733133d-26, 8.81804906d-27,&
      7.00417973d-27, 5.43356567d-27, 4.09857884d-27, 2.99651764d-27,&
      2.11902962d-27, 1.45014127d-27, 9.62291023d-28, 6.21548647d-28,&
      3.92807578d-28, 2.44230375d-28, 1.50167782d-28, 9.17611405d-29,&
      5.58707641d-29, 3.40570915d-29, 2.08030862d-29, 1.27588676d-29,&
      7.86535588d-30, 4.87646151d-30, 3.03888897d-30, 1.90578649d-30,&
      1.20195947d-30, 7.61955060d-31, 4.85602199d-31, 3.11049969d-31,&
      2.00087065d-31, 1.29223740d-31, 8.37422008d-32, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00 /

  data f_255 /  0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 1.76014287d-44,&
      5.07057938d-42, 7.90473970d-40, 7.31852999d-38, 4.30709255d-36,&
      1.70009061d-34, 4.67925160d-33, 9.21703546d-32, 1.31918676d-30,&
      1.37393161d-29, 1.03102379d-28, 5.45694018d-28, 1.98699648d-27,&
      4.90346776d-27, 8.40524725d-27, 1.08321456d-26, 1.15714525d-26,&
      1.10905152d-26, 1.00155023d-26, 8.75799161d-27, 7.50935839d-27,&
      6.35253533d-27, 5.30919268d-27, 4.37669455d-27, 3.55185164d-27,&
      2.82347055d-27, 2.19257595d-27, 1.65589541d-27, 1.21224987d-27,&
      8.58395132d-28, 5.88163935d-28, 3.90721447d-28, 2.52593407d-28,&
      1.59739995d-28, 9.93802874d-29, 6.11343388d-29, 3.73711135d-29,&
      2.27618743d-29, 1.38793199d-29, 8.48060787d-30, 5.20305940d-30,&
      3.20867365d-30, 1.99011277d-30, 1.24064551d-30, 7.78310544d-31,&
      4.91013681d-31, 3.11338381d-31, 1.98451675d-31, 1.27135460d-31,&
      8.17917486d-32, 5.28280497d-32, 3.42357159d-32, 0.00000000d+00,&
      0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00 /

  abstract interface
    subroutine get_subr1(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,res)
      use mod_global_parameters
      integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in) :: w(ixImin1:ixImax1,nw)
      double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(out):: res(ixImin1:ixImax1)
    end subroutine get_subr1

  end interface

  type te_fluid

    procedure (get_subr1), pointer, nopass :: get_rho => null()
    procedure (get_subr1), pointer, nopass :: get_pthermal => null()
    procedure (get_subr1), pointer, nopass :: get_var_Rfactor => null()

    ! factor in eq of state p = Rfactor * rho * T
    ! used for getting temperature
    double precision :: Rfactor = 1d0

  end type te_fluid


  contains

    subroutine get_line_info(wl,ion,mass,logTe,line_center,spatial_px,&
       spectral_px,sigma_PSF,width_slit)
      ! get information of the spectral line
      ! wl: wavelength
      ! mass: ion mass, unit -- proton mass
      ! logTe: peak temperature of emission line in logarithm
      ! line_center: center wavelength of emission line, unit -- Angstrom (0.1 nm) 
      ! spatial_px: pixel size in space of instrument (for image), unit -- arcsec
      ! spectral_px: pixel size in wagelength of instrument (for spectrum), unit -- Angstrom
      ! sigma_PSF: width of point spread function core (for instrument), unit -- pixel
      ! width_slit: width of slit for spectrograph, unit -- arcsec
      use mod_global_parameters

      integer, intent(in) :: wl
      integer, intent(out) :: mass
      character(len=30), intent(out) :: ion
      double precision, intent(out) :: logTe,line_center,spatial_px,&
         spectral_px
      double precision, intent(out) :: sigma_PSF,width_slit

      select case(wl)
      case(304)
        ion='He II'
        mass=4
        logTe=4.7d0
        line_center=303.8d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=0.895d0
        width_slit=0.6d0
      case(171)
        ion='Fe IX'
        mass=56
        logTe=5.8d0
        line_center=171.1d0
        spatial_px=0.6d0
        spectral_px=0.02d0 
        sigma_PSF=1.019d0
        width_slit=0.6d0
      case(193)
        ion='Fe XXIV'
        mass=56
        logTe=7.3d0
        line_center=193.5d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=0.813d0
        width_slit=0.6d0
      case(211)
        ion='Fe XIV'
        mass=56
        logTe=6.3d0
        line_center=211.3d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=0.913d0
        width_slit=0.6d0
      case(335)
        ion='Fe XVI'
        mass=56
        logTe=6.4d0
        line_center=335.4d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=1.019d0
        width_slit=0.6d0
      case(94)
        ion='Fe XVIII'
        mass=56
        logTe=6.8d0
        line_center=93.9d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=1.025d0
        width_slit=0.6d0
      case(131)
        ion='Fe XXI'
        mass=56
        logTe=7.0d0
        line_center=131.0d0
        spatial_px=0.6d0
        spectral_px=0.02d0
        sigma_PSF=0.984d0
        width_slit=0.6d0
      case(1354)
        ion='Fe XXI'
        mass=56
        logTe=7.0d0
        line_center=1354.1d0
        spatial_px=0.1663d0
        spectral_px=12.98d-3
        sigma_PSF=1.d0
        width_slit=0.33d0
      case(263)
        ion='Fe XVI'
        mass=56
        logTe=6.4d0
        line_center=262.976d0
        spatial_px=1.d0
        spectral_px=22.d-3
        sigma_PSF=1.d0
        width_slit=2.d0
      case(264)
        ion='Fe XXIII'
        mass=56
        logTe=7.1d0
        line_center=263.765d0
        spatial_px=1.d0
        spectral_px=22.d-3
        sigma_PSF=1.d0
        width_slit=2.d0
      case(192)
        ion='Fe XXIV'
        mass=56
        logTe=7.2d0
        line_center=192.028d0
        spatial_px=1.d0
        spectral_px=22.d-3
        sigma_PSF=1.d0
        width_slit=2.d0
      case(255)
        ion='Fe XXIV'
        mass=56
        logTe=7.2d0
        line_center=255.113d0
        spatial_px=1.d0
        spectral_px=22.d-3
        sigma_PSF=1.d0
        width_slit=2.d0
      case default
        call mpistop("No information about this line")
      end select
    end subroutine get_line_info
    
    subroutine get_EUV(wl,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,fl,flux)
      ! calculate the local emission intensity of given EUV line (optically thin)
      ! wavelength is the wave length of the emission line
      ! unit [DN cm^-1 s^-1 pixel^-1]
      ! ingrate flux along line of sight: DN s^-1 pixel^-1
      use mod_global_parameters

      integer, intent(in) :: wl
      integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
      double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in) :: w(ixImin1:ixImax1,1:nw)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: flux(ixImin1:ixImax1)

      integer :: n_table
      double precision, allocatable :: t_table(:),f_table(:)
      integer :: ix1,iTt,i
      double precision :: pth(ixImin1:ixImax1),Te(ixImin1:ixImax1),&
         Ne(ixImin1:ixImax1)
      double precision :: logT,logGT,GT

      ! selecting emission table 
      select case(wl)
      case(94)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_94(1:n_aia)
      case(131)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_131(1:n_aia)
      case(171)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_171(1:n_aia)
      case(193)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_193(1:n_aia)
      case(211)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_211(1:n_aia)
      case(304)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_304(1:n_aia)
      case(335)
        n_table=n_aia
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_aia(1:n_aia)
        f_table(1:n_table)=f_335(1:n_aia)
      case(1354)
        n_table=n_iris
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_iris(1:n_iris)
        f_table(1:n_table)=f_1354(1:n_iris)
      case(263)
        n_table=n_eis
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_eis1(1:n_eis)
        f_table(1:n_table)=f_263(1:n_eis)
      case(264)
        n_table=n_eis
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_eis2(1:n_eis)
        f_table(1:n_table)=f_264(1:n_eis)
      case(192)
        n_table=n_eis
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_eis2(1:n_eis)
        f_table(1:n_table)=f_192(1:n_eis)
      case(255)
        n_table=n_eis
        allocate(t_table(1:n_table))
        allocate(f_table(1:n_table))
        t_table(1:n_table)=t_eis2(1:n_eis)
        f_table(1:n_table)=f_255(1:n_eis)
      case default
        allocate(t_table(1))
        allocate(f_table(1))
        call mpistop("Unknown wavelength")
      end select
      call fl%get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)
      call fl%get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Ne)
      if(associated(fl%get_var_Rfactor)) then
        call fl%get_var_Rfactor(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Te)
        Te(ixOmin1:ixOmax1)=pth(ixOmin1:ixOmax1)/(Ne(ixOmin1:ixOmax1)*Te(&
           ixOmin1:ixOmax1))*unit_temperature
      else  
        Te(ixOmin1:ixOmax1)=pth(ixOmin1:ixOmax1)/(Ne(&
           ixOmin1:ixOmax1)*fl%Rfactor)*unit_temperature
      endif
      if (SI_unit) then
        Ne(ixOmin1:ixOmax1)=Ne(ixOmin1:ixOmax1)*unit_numberdensity/1.d6 !m-3 -> cm-3
        flux(ixOmin1:ixOmax1)=Ne(ixOmin1:ixOmax1)**2
      else
        Ne(ixOmin1:ixOmax1)=Ne(ixOmin1:ixOmax1)*unit_numberdensity
        flux(ixOmin1:ixOmax1)=Ne(ixOmin1:ixOmax1)**2
      endif

      select case(wl)
      case(94,131,171,193,211,304,335,1354)
      ! temperature table in log
        do i=1,n_table
          if(f_table(i) .gt. 1.d-99) then
            f_table(i)=log10(f_table(i))
          else
            f_table(i)=-99.d0
          endif
        enddo 
        logGT=zero
        do ix1=ixOmin1,ixOmax1
          logT=log10(Te(ix1))
          if (logT>=t_table(1) .and. logT<=t_table(n_table)) then
            do iTt=1,n_table-1
              if (logT>=t_table(iTt) .and. logT<t_table(iTt+1)) then
                logGT=f_table(iTt)*(logT-t_table(iTt+&
                   1))/(t_table(iTt)-t_table(iTt+1))+f_table(iTt+&
                   1)*(logT-t_table(iTt))/(t_table(iTt+1)-t_table(iTt))
              endif
            enddo
            flux(ix1)=flux(ix1)*(10**(logGT))
            if(flux(ix1)<smalldouble) flux(ix1)=0.d0
          else
            flux(ix1)=zero
          endif
        enddo
      case default
      ! temperature table linear
        do ix1=ixOmin1,ixOmax1
          if (Te(ix1)>=t_table(1) .and. Te(ix1)<=t_table(n_table)) then
            do iTt=1,n_table-1
              if (Te(ix1)>=t_table(iTt) .and. Te(ix1)<t_table(iTt+1)) then
                GT=f_table(iTt)*(Te(ix1)-t_table(iTt+&
                   1))/(t_table(iTt)-t_table(iTt+1))+f_table(iTt+&
                   1)*(Te(ix1)-t_table(iTt))/(t_table(iTt+1)-t_table(iTt))
              endif
            enddo
            flux(ix1)=flux(ix1)*GT
            if(flux(ix1)<smalldouble) flux(ix1)=0.d0
          else
            flux(ix1)=zero
          endif
        enddo
      end select

      deallocate(t_table,f_table)
    end subroutine get_EUV

    subroutine get_SXR(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,fl,flux,El,Eu)
      !synthesize thermal SXR from El keV to Eu keV released by cm^-3/m^-3
      ! volume of plasma in 1 s
      !flux (cgs): photons cm^-3 s^-1
      !flux (SI): photons m^-3 s^-1
      use mod_global_parameters

      integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
      integer, intent(in)           :: El,Eu
      double precision, intent(in)  :: x(ixImin1:ixImax1,1:ndim)
      double precision, intent(in)  :: w(ixImin1:ixImax1,nw)
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: flux(ixImin1:ixImax1)

      integer :: ix1,ixO1
      integer :: iE,numE
      double precision :: I0,kb,keV,dE,Ei
      double precision :: pth(ixImin1:ixImax1),Te(ixImin1:ixImax1),&
         kbT(ixImin1:ixImax1)
      double precision :: Ne(ixImin1:ixImax1),gff(ixImin1:ixImax1),&
         fi(ixImin1:ixImax1)
      double precision :: EM(ixImin1:ixImax1)

      I0=3.01d-15   ! I0*4*pi*AU**2, I0 from Pinto (2015)
      kb=const_kb
      keV=1.0d3*const_ev
      dE=0.1
      numE=floor((Eu-El)/dE)
      call fl%get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)
      call fl%get_rho(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Ne)

      if(associated(fl%get_var_Rfactor)) then
        call fl%get_var_Rfactor(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Te)
        Te(ixOmin1:ixOmax1)=pth(ixOmin1:ixOmax1)/(Ne(ixOmin1:ixOmax1)*Te(&
           ixOmin1:ixOmax1))*unit_temperature
      else
        Te(ixOmin1:ixOmax1)=pth(ixOmin1:ixOmax1)/(Ne(&
           ixOmin1:ixOmax1)*fl%Rfactor)*unit_temperature
      endif
      if (SI_unit) then
        Ne(ixOmin1:ixOmax1)=Ne(ixOmin1:ixOmax1)*unit_numberdensity/1.d6 !m-3 -> cm-3
        EM(ixOmin1:ixOmax1)=(Ne(ixOmin1:ixOmax1))**2*1.d6 ! cm-3 m-3
      else
        Ne(ixOmin1:ixOmax1)=Ne(ixOmin1:ixOmax1)*unit_numberdensity
        EM(ixOmin1:ixOmax1)=(Ne(ixOmin1:ixOmax1))**2
      endif
      kbT(ixOmin1:ixOmax1)=kb*Te(ixOmin1:ixOmax1)/keV
      flux(ixOmin1:ixOmax1)=0.0d0
      do iE=0,numE-1
        Ei=dE*iE+El*1.d0
        gff(ixOmin1:ixOmax1)=1.d0
        do ix1=ixOmin1,ixOmax1
          if (kbT(ix1)>0.01*Ei) then
            if(kbT(ix1)<Ei) gff(ix1)=(kbT(ix1)/Ei)**0.4
            fi(ix1)=(EM(ix1)*gff(ix1))*dexp(-&
               Ei/(kbT(ix1)))/(Ei*dsqrt(kbT(ix1)))
          else
            fi(ix1)=zero
          endif
        enddo
        flux(ixOmin1:ixOmax1)=flux(ixOmin1:ixOmax1)+fi(ixOmin1:ixOmax1)*dE
      enddo
      flux(ixOmin1:ixOmax1)=flux(ixOmin1:ixOmax1)*I0
    end subroutine get_SXR

    subroutine get_GOES_SXR_flux(xboxmin1,xboxmax1,fl,eflux)
      !get GOES SXR 1-8A flux observing at 1AU from given box [w/m^2]
      use mod_global_parameters

      double precision, intent(in) :: xboxmin1,xboxmax1
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: eflux

      double precision :: dxb1,xbmin1,xbmax1
      integer :: iigrid,igrid,j
      integer :: ixOmin1,ixOmax1,ixImin1,ixImax1,ix1
      double precision :: eflux_grid,eflux_pe

      ixImin1=ixglo1;
      ixImax1=ixghi1;
      ixOmin1=ixmlo1;
      ixOmax1=ixmhi1;
      eflux_pe=zero
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        dxlevel(1)=rnode(rpdx1_,igrid);
        xbmin1=rnode(rpxmin1_,igrid);
        xbmax1=rnode(rpxmax1_,igrid);
        call get_GOES_flux_grid(ixImin1,ixImax1,ixOmin1,ixOmax1,ps(igrid)%w,&
           ps(igrid)%x,ps(igrid)%dvolume(ixImin1:ixImax1),xboxmin1,xboxmax1,&
           xbmin1,xbmax1,fl,eflux_grid)
        eflux_pe=eflux_pe+eflux_grid
      enddo
      call MPI_ALLREDUCE(eflux_pe,eflux,1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,&
         ierrmpi)
    end subroutine get_GOES_SXR_flux

    subroutine get_GOES_flux_grid(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,dV,&
       xboxmin1,xboxmax1,xbmin1,xbmax1,fl,eflux_grid)
      use mod_global_parameters

      integer, intent(in)           :: ixImin1,ixImax1,ixOmin1,ixOmax1
      double precision, intent(in)  :: x(ixImin1:ixImax1,1:ndim),&
         dV(ixImin1:ixImax1)
      double precision, intent(in)  :: w(ixImin1:ixImax1,nw)
      double precision, intent(in)  :: xboxmin1,xboxmax1,xbmin1,xbmax1
      type(te_fluid), intent(in) :: fl
      double precision, intent(out) :: eflux_grid

      integer :: ix1,ixO1,ixbmin1,ixbmax1
      integer :: iE,numE,j,inbox
      double precision :: I0,kb,keV,dE,Ei,El,Eu,A_cgs
      double precision :: pth(ixImin1:ixImax1),Te(ixImin1:ixImax1),&
         kbT(ixImin1:ixImax1)
      double precision :: Ne(ixImin1:ixImax1),EM(ixImin1:ixImax1)
      double precision :: gff,fi,erg_SI

      ! check whether the grid is inside given box
      inbox=0
      if (xbmin1<xboxmax1 .and. xbmax1>xboxmin1) inbox=inbox+1

      if (inbox==ndim) then
        ! indexes for cells inside given box
        ixbmin1=ixOmin1;
        ixbmax1=ixOmax1;
        if (xbmax1>xboxmax1) ixbmax1=ixOmax1-ceiling((xbmax1-&
           xboxmax1)/dxlevel(1))
        if (xbmin1<xboxmin1) ixbmin1=ceiling((xboxmin1-xbmin1)/dxlevel(1))+&
           ixOmin1

        I0=1.07d-38 !photon flux index for observed at 1AU [photon cm3 m-2 s-1 keV-1]
        kb=const_kb
        keV=1.0d3*const_ev
        erg_SI=1.d-7
        A_cgs=1.d-8 ! Angstrom
        El=const_h*const_c/(8.d0*A_cgs)/keV ! 8 A
        Eu=const_h*const_c/(1.d0*A_cgs)/keV ! 1 A
        dE=0.1  ! keV
        numE=floor((Eu-El)/dE)
        call fl%get_pthermal(w,x,ixImin1,ixImax1,ixbmin1,ixbmax1,pth)
        call fl%get_rho(w,x,ixImin1,ixImax1,ixbmin1,ixbmax1,Ne)
        if(associated(fl%get_var_Rfactor)) then
          call fl%get_var_Rfactor(w,x,ixImin1,ixImax1,ixbmin1,ixbmax1,Te)
          Te(ixbmin1:ixbmax1)=pth(ixbmin1:ixbmax1)/(Ne(ixbmin1:ixbmax1)*Te(&
             ixbmin1:ixbmax1))*unit_temperature
        else
          Te(ixbmin1:ixbmax1)=pth(ixbmin1:ixbmax1)/(Ne(&
             ixbmin1:ixbmax1)*fl%Rfactor)*unit_temperature
        endif
        if (SI_unit) then
          Ne(ixOmin1:ixOmax1)=Ne(ixOmin1:ixOmax1)*unit_numberdensity/1.d6 !m-3 -> cm-3
          EM(ixbmin1:ixbmax1)=(I0*(Ne(ixbmin1:ixbmax1))**&
             2)*dV(ixbmin1:ixbmax1)*(unit_length*1.d2)**3 !cm-3
        else
          Ne(ixOmin1:ixOmax1)=Ne(ixOmin1:ixOmax1)*unit_numberdensity
          EM(ixbmin1:ixbmax1)=(I0*(Ne(ixbmin1:ixbmax1))**&
             2)*dV(ixbmin1:ixbmax1)*unit_length**3
        endif
        kbT(ixbmin1:ixbmax1)=kb*Te(ixbmin1:ixbmax1)/keV
        eflux_grid=0.0d0

        do iE=0,numE-1
          Ei=dE*iE+El
          do ix1=ixbmin1,ixbmax1
            if (kbT(ix1)>1.d-2*Ei) then
              if(kbT(ix1)<Ei) then
                gff=(kbT(ix1)/Ei)**0.4
              else
                gff=1.d0
              endif
              fi=(EM(ix1)*gff)*dexp(-Ei/(kbT(ix1)))/(Ei*dsqrt(kbT(ix1)))
              eflux_grid=eflux_grid+fi*dE*Ei
            endif
          enddo
        enddo
        eflux_grid=eflux_grid*keV*erg_SI
      endif

    end subroutine get_GOES_flux_grid

  

    subroutine write_image_vtiCC(qunit,xO1,xO2,dxO1,dxO2,wO,nXO1,nXO2,nWO,nC1,&
       nC2)
      ! write image data to vti
      use mod_global_parameters

      integer, intent(in) :: qunit,nXO1,nXO2,nWO,nC1,nC2
      double precision, intent(in) :: xO1(nXO1),xO2(nxO2)
      double precision, intent(in) :: dxO1(nxO1),dxO2(nxO2)
      double precision, intent(in) :: wO(nXO1,nXO2,nWO)

      double precision :: origin(1:3), spacing(1:3)
      integer :: wholeExtent(1:6),extent(1:6)
      integer :: nP1,nP2,iP1,iP2,iw
      integer :: ixC1,ixC2,ixCmin1,ixCmax1,ixCmin2,ixCmax2

      integer :: filenr
      logical :: fileopen
      character (70) :: subname,wname,vname,nameL,nameS
      character (len=std_len) :: filename
      integer :: mass
      double precision :: logTe
      character (30) :: ion
      double precision :: line_center
      double precision :: spatial_rsl,spectral_rsl,sigma_PSF,wslit


      origin(1)=xO1(1)-0.5d0*dxO1(1)
      origin(2)=xO2(1)-0.5d0*dxO2(1)
      origin(3)=zero
      spacing(1)=dxO1(1)
      spacing(2)=dxO2(1)
      spacing(3)=zero
      wholeExtent=zero
      wholeExtent(2)=nXO1
      wholeExtent(4)=nXO2
      nP1=nXO1/nC1
      nP2=nXO2/nC2

      ! get information of emission line
      if (convert_type=='EIvtiCCmpi') then
        call get_line_info(wavelength,ion,mass,logTe,line_center,spatial_rsl,&
           spectral_rsl,sigma_PSF,wslit)
      else if (convert_type=='ESvtiCCmpi') then
        call get_line_info(spectrum_wl,ion,mass,logTe,line_center,spatial_rsl,&
           spectral_rsl,sigma_PSF,wslit)
      endif

      if (mype==0) then
        inquire(qunit,opened=fileopen)
        if(.not.fileopen)then
          ! generate filename 
          filenr=snapshotini
          if (autoconvert) filenr=snapshotnext
          if (convert_type=='EIvtiCCmpi') then
            write(filename,'(a,i4.4,a)') trim(filename_euv),filenr,".vti"
          else if (convert_type=='SIvtiCCmpi') then
            write(filename,'(a,i4.4,a)') trim(filename_sxr),filenr,".vti"
          else if (convert_type=='ESvtiCCmpi') then
            write(filename,'(a,i4.4,a)') trim(filename_spectrum),filenr,".vti"
          endif
          open(qunit,file=filename,status='unknown',form='formatted')
        endif

        ! generate xml header
        write(qunit,'(a)')'<?xml version="1.0"?>'
        write(qunit,'(a)',advance='no') '<VTKFile type="ImageData"'
        write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
        write(qunit,'(a,3(1pe14.6),a,6(i10),a,3(1pe14.6),a)')&
           '  <ImageData Origin="',origin,'" WholeExtent="',wholeExtent,&
           '" Spacing="',spacing,'">'
        ! file info        
        write(qunit,'(a)')'<FieldData>'
        write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
           'NumberOfTuples="1" format="ascii">'
        write(qunit,*) real(global_time*time_convert_factor)
        write(qunit,'(a)')'</DataArray>'
        if (convert_type=='EIvtiCCmpi' .or. convert_type=='ESvtiCCmpi') then
          write(qunit,'(2a)')'<DataArray type="Float32" Name="logT" ',&
             'NumberOfTuples="1" format="ascii">'
          write(qunit,*) real(logTe)
          write(qunit,'(a)')'</DataArray>'
        endif
        write(qunit,'(a)')'</FieldData>'
        ! pixel/cell data
        do iP1=1,nP1
          do iP2=1,nP2
            extent=zero
            extent(1)=(iP1-1)*nC1
            extent(2)=iP1*nC1
            extent(3)=(iP2-1)*nC2
            extent(4)=iP2*nC2
            ixCmin1=extent(1)+1
            ixCmax1=extent(2)
            ixCmin2=extent(3)+1
            ixCmax2=extent(4)
            write(qunit,'(a,6(i10),a)') '<Piece Extent="',extent,'">'
            write(qunit,'(a)')'<CellData>'
            do iw=1,nWO
              ! variable name
              if (convert_type=='EIvtiCCmpi') then
                if (iw==1) then
                  if (wavelength<100) then
                    write(vname,'(a,i2)') "AIA",wavelength
                  else if (wavelength<1000) then
                    write(vname,'(a,i3)') "AIA",wavelength
                  else
                    write(vname,'(a,i4)') "IRIS",wavelength
                  endif
                endif
                if (iw==2) vname='Doppler_velocity'
              else if (convert_type=='SIvtiCCmpi') then
                if (emin_sxr<10 .and. emax_sxr<10) then
                  write(vname,'(a,i1,a,i1,a)') "SXR",emin_sxr,"-",emax_sxr,&
                     "keV"
                else if (emin_sxr<10 .and. emax_sxr>=10) then
                  write(vname,'(a,i1,a,i2,a)') "SXR",emin_sxr,"-",emax_sxr,&
                     "keV"
                else
                  write(vname,'(a,i2,a,i2,a)') "SXR",emin_sxr,"-",emax_sxr,&
                     "keV"
                endif
              else if (convert_type=='ESvtiCCmpi') then
                if (spectrum_wl==1354) then
                  write(vname,'(a,i4)') "SG",spectrum_wl
                else
                  write(vname,'(a,i3)') "EIS",spectrum_wl
                endif
              endif
              write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
                 TRIM(vname),'" format="ascii">'
              write(qunit,'(200(1pe14.6))') ((wO(ixC1,ixC2,iw),ixC1=ixCmin1,&
                 ixCmax1),ixC2=ixCmin2,ixCmax2)
              write(qunit,'(a)')'</DataArray>'
            enddo
            write(qunit,'(a)')'</CellData>'
            write(qunit,'(a)')'</Piece>'
          enddo
        enddo
        ! end
        write(qunit,'(a)')'</ImageData>'
        write(qunit,'(a)')'</VTKFile>'
        close(qunit)
      endif

    end subroutine write_image_vtiCC

    subroutine write_image_vtuCC(qunit,xC,wC,dxC,nPiece,nC1,nC2,nWC,datatype)
      ! write image data to vtu
      use mod_global_parameters

      integer, intent(in) :: qunit
      integer, intent(in) :: nPiece,nC1,nC2,nWC
      double precision, intent(in) :: xC(nPiece,nC1,nC2,2),dxC(nPiece,nc1,nc2,&
         2)
      double precision, intent(in) :: wC(nPiece,nC1,nC2,nWC)
      character(20), intent(in) :: datatype

      integer :: nP1,nP2
      double precision :: xP(nPiece,nC1+1,nC2+1,2)
      integer :: filenr
      logical :: fileopen
      character (70) :: subname,wname,vname,nameL,nameS
      character (len=std_len) :: filename
      integer :: ixC1,ixC2,ixP,ix1,ix2,j
      integer :: nc,np,icel,VTK_type
      integer :: mass
      double precision :: logTe
      character (30) :: ion
      double precision :: line_center
      double precision :: spatial_rsl,spectral_rsl,sigma_PSF,wslit

      nP1=nC1+1
      nP2=nC2+1
      np=nP1*nP2
      nc=nC1*nC2
      ! cell corner location     
      do ixP=1,nPiece
        do ix1=1,nP1
          do ix2=1,nP2
            if (ix1<nP1) xP(ixP,ix1,ix2,1)=xC(ixP,ix1,1,1)-0.5d0*dxC(ixP,ix1,1,&
               1)
            if (ix1==nP1) xP(ixP,ix1,ix2,1)=xC(ixP,ix1-1,1,1)+0.5d0*dxC(ixP,&
               ix1-1,1,1)
            if (ix2<nP2) xP(ixP,ix1,ix2,2)=xC(ixP,1,ix2,2)-0.5d0*dxC(ixP,1,ix2,&
               2)
            if (ix2==nP2) xP(ixP,ix1,ix2,2)=xC(ixP,1,ix2-1,2)+0.5d0*dxC(ixP,1,&
               ix2-1,2)
          enddo
        enddo
      enddo
      ! get information of emission line
      if (datatype=='image_euv') then
        call get_line_info(wavelength,ion,mass,logTe,line_center,spatial_rsl,&
           spectral_rsl,sigma_PSF,wslit)
      else if (datatype=='spectrum_euv') then
        call get_line_info(spectrum_wl,ion,mass,logTe,line_center,spatial_rsl,&
           spectral_rsl,sigma_PSF,wslit)
      endif

      if (mype==0) then
        inquire(qunit,opened=fileopen)
        if(.not.fileopen)then
          ! generate filename 
          filenr=snapshotini
          if (autoconvert) filenr=snapshotnext
          if (datatype=='image_euv') then
            write(filename,'(a,i4.4,a)') trim(filename_euv),filenr,".vtu"
          else if (datatype=='image_sxr') then
            write(filename,'(a,i4.4,a)') trim(filename_sxr),filenr,".vtu"
          else if (datatype=='spectrum_euv') then
            write(filename,'(a,i4.4,a)') trim(filename_spectrum),filenr,".vtu"
          endif
          open(qunit,file=filename,status='unknown',form='formatted')
        endif
        ! generate xml header
        write(qunit,'(a)')'<?xml version="1.0"?>'
        write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
        write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
        write(qunit,'(a)')'<UnstructuredGrid>'
        write(qunit,'(a)')'<FieldData>'
        write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
           'NumberOfTuples="1" format="ascii">'
        write(qunit,*) real(global_time*time_convert_factor)
        write(qunit,'(a)')'</DataArray>'
        if (datatype=='image_euv' .or. datatype=='spectrum_euv') then
          write(qunit,'(2a)')'<DataArray type="Float32" Name="logT" ',&
             'NumberOfTuples="1" format="ascii">'
          write(qunit,*) real(logTe)
          write(qunit,'(a)')'</DataArray>'
        endif
        write(qunit,'(a)')'</FieldData>'
        do ixP=1,nPiece
          write(qunit,'(a,i7,a,i7,a)') '<Piece NumberOfPoints="',np,&
             '" NumberOfCells="',nc,'">'
          write(qunit,'(a)')'<CellData>'
          do j=1,nWC
            if (datatype=='image_euv') then
              if (j==1) then
                if (wavelength<100) then
                  write(vname,'(a,i2)') "AIA",wavelength
                else if (wavelength<1000) then
                  write(vname,'(a,i3)') "AIA",wavelength
                else
                  write(vname,'(a,i4)') "IRIS",wavelength
                endif
              endif
              if (j==2) vname='Doppler_velocity'
            else if (datatype=='image_sxr') then
              if (emin_sxr<10 .and. emax_sxr<10) then
                write(vname,'(a,i1,a,i1,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
              else if (emin_sxr<10 .and. emax_sxr>=10) then
                write(vname,'(a,i1,a,i2,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
              else
                write(vname,'(a,i2,a,i2,a)') "SXR",emin_sxr,"-",emax_sxr,"keV"
              endif
            else if (datatype=='spectrum_euv') then
              if (spectrum_wl==1354) then
                write(vname,'(a,i4)') "SG",spectrum_wl
              else
                write(vname,'(a,i3)') "EIS",spectrum_wl
              endif
            endif
            write(qunit,'(a,a,a)')'<DataArray type="Float64" Name="',&
               TRIM(vname),'" format="ascii">'
            write(qunit,'(200(1pe14.6))') ((wC(ixP,ixC1,ixC2,j),ixC1=1,nC1),&
               ixC2=1,nC2)
            write(qunit,'(a)')'</DataArray>'
          enddo
          write(qunit,'(a)')'</CellData>'
          write(qunit,'(a)')'<Points>'
          write(qunit,'(a)'&
)'<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          do ix2=1,nP2
            do ix1=1,nP1 
              if (datatype=='image_euv' .and. resolution_euv=='data') then
                if (LOS_phi==0 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') 0.d0,xP(ixP,ix1,ix2,1),xP(ixP,&
                     ix1,ix2,2)
                else if (LOS_phi==90 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,2),0.d0,xP(ixP,&
                     ix1,ix2,1)
                else
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,&
                     2),0.d0
                endif
              else if (datatype=='image_sxr' .and. resolution_sxr=='data') &
                 then
                if (LOS_phi==0 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') 0.d0,xP(ixP,ix1,ix2,1),xP(ixP,&
                     ix1,ix2,2)
                else if (LOS_phi==90 .and. LOS_theta==90) then
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,2),0.d0,xP(ixP,&
                     ix1,ix2,1)
                else
                  write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,&
                     2),0.d0
                endif
              else
                write(qunit,'(3(1pe14.6))') xP(ixP,ix1,ix2,1),xP(ixP,ix1,ix2,&
                   2),0.d0
              endif            
            enddo
          enddo
          write(qunit,'(a)')'</DataArray>'
          write(qunit,'(a)')'</Points>'
          ! connetivity part
          write(qunit,'(a)')'<Cells>'
          write(qunit,'(a)'&
             )'<DataArray type="Int32" Name="connectivity" format="ascii">'
          do ix2=1,nC2
            do ix1=1,nC1
              write(qunit,'(4(i7))') ix1-1+(ix2-1)*nP1,ix1+(ix2-1)*nP1,&
                 ix1-1+ix2*nP1,ix1+ix2*nP1
            enddo
          enddo
          write(qunit,'(a)')'</DataArray>'
          ! offsets data array
          write(qunit,'(a)'&
             )'<DataArray type="Int32" Name="offsets" format="ascii">'
          do icel=1,nc
            write(qunit,'(i7)') icel*(2**2)
          enddo
          write(qunit,'(a)')'</DataArray>'
          ! VTK cell type data array
          write(qunit,'(a)'&
             )'<DataArray type="Int32" Name="types" format="ascii">'
          ! VTK_LINE=3; VTK_PIXEL=8; VTK_VOXEL=11 -> vtk-syntax
          VTK_type=8        
          do icel=1,nc
            write(qunit,'(i2)') VTK_type
          enddo
          write(qunit,'(a)')'</DataArray>' 
          write(qunit,'(a)')'</Cells>'
          write(qunit,'(a)')'</Piece>'
        enddo
        write(qunit,'(a)')'</UnstructuredGrid>'
        write(qunit,'(a)')'</VTKFile>'
        close(qunit)
      endif
    end subroutine write_image_vtuCC

    subroutine dot_product_loc(vec_in1,vec_in2,res)
      double precision, intent(in) :: vec_in1(1:3),vec_in2(1:3)
      double precision, intent(out) :: res
      integer :: j

      res=zero
      do j=1,3
        res=res+vec_in1(j)*vec_in2(j)
      enddo

    end subroutine dot_product_loc

    subroutine cross_product_loc(vec_in1,vec_in2,vec_out)
      double precision, intent(in) :: vec_in1(1:3),vec_in2(1:3)
      double precision, intent(out) :: vec_out(1:3)

      vec_out(1)=vec_in1(2)*vec_in2(3)-vec_in1(3)*vec_in2(2)
      vec_out(2)=vec_in1(3)*vec_in2(1)-vec_in1(1)*vec_in2(3)
      vec_out(3)=vec_in1(1)*vec_in2(2)-vec_in1(2)*vec_in2(1)

    end subroutine cross_product_loc

    subroutine init_vectors()
      integer :: j
      double precision :: LOS_psi
      double precision :: vec_z(1:3),vec_temp1(1:3),vec_temp2(1:3)

      ! vectors for image coordinate
      vec_LOS(1)=-cos(dpi*LOS_phi/180.d0)*sin(dpi*LOS_theta/180.d0)
      vec_LOS(2)=-sin(dpi*LOS_phi/180.d0)*sin(dpi*LOS_theta/180.d0)
      vec_LOS(3)=-cos(dpi*LOS_theta/180.d0)
      do j=1,3
        if (abs(vec_LOS(j))<=smalldouble) vec_LOS(j)=zero
      enddo
      vec_z(:)=zero
      vec_z(3)=1.d0
      if (LOS_theta==zero) then
        vec_xI1=zero
        vec_xI2=zero
        vec_xI1(1)=1.d0
        vec_xI2(2)=1.d0
      else
        call cross_product_loc(vec_LOS,vec_z,vec_xI1)
        call cross_product_loc(vec_xI1,vec_LOS,vec_xI2)
      endif
      vec_temp1=vec_xI1/sqrt(vec_xI1(1)**2+vec_xI1(2)**2+vec_xI1(3)**2)
      vec_temp2=vec_xI2/sqrt(vec_xI2(1)**2+vec_xI2(2)**2+vec_xI2(3)**2)
      LOS_psi=dpi*image_rotate/180.d0
      vec_xI1=vec_temp1*cos(LOS_psi)-vec_temp2*sin(LOS_psi)
      vec_xI2=vec_temp2*cos(LOS_psi)+vec_temp1*sin(LOS_psi)

      do j=1,3
        if (abs(vec_xI1(j))<smalldouble) vec_xI1(j)=zero
        if (abs(vec_xI2(j))<smalldouble) vec_xI2(j)=zero
      enddo

      if (mype==0) write(*,'(a,f5.2,f6.2,f6.2,a)') ' LOS vector: [',vec_LOS(1),&
         vec_LOS(2),vec_LOS(3),']'
      if (mype==0) write(*,'(a,f5.2,f6.2,f6.2,a)') ' xI1 vector: [',vec_xI1(1),&
         vec_xI1(2),vec_xI1(3),']'
      if (mype==0) write(*,'(a,f5.2,f6.2,f6.2,a)') ' xI2 vector: [',vec_xI2(1),&
         vec_xI2(2),vec_xI2(3),']'

    end subroutine init_vectors

    subroutine get_cor_image(x_3D,x_image)
      double precision :: x_3D(1:3),x_image(1:2)
      double precision :: res,res_origin

      call dot_product_loc(x_3D,vec_xI1,res)
      call dot_product_loc(x_origin,vec_xI1,res_origin)
      x_image(1)=res-res_origin
      call dot_product_loc(x_3D,vec_xI2,res)
      call dot_product_loc(x_origin,vec_xI2,res_origin)
      x_image(2)=res-res_origin

    end subroutine get_cor_image

end module mod_thermal_emission
