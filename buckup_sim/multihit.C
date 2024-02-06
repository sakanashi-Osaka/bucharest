{

  //  tree->Draw("mK3a[]+mK3a[(Iteration$+1)%3]:mK3a[]+mK3a[(Iteration$+1)%3]+mK4He>>h(100,14,18,100,0,5)","m3achF[(Iteration$+0)%3]==m3achF[(Iteration$+1)%3] && m3achR[(Iteration$+0)%3]==m3achR[(Iteration$+1)%3] && idp4He && idp3a[(Iteration$+0)%3] && idp3a[(Iteration$+1)%3] && m3aseg[(Iteration$+2)%3]!=((m4Heseg+2)%4) && abs(abs(mph4He-mph3a[])-180)<30 && mth4He>-2*mth3a[]+120","col"); //alpha 2hit

  tree->Draw("mK12C:mK12C+mK4He>>hh(100,14,18,100,0,5)"," idp4He && idp12C && abs(abs(mph4He-mph3a[])-180)<30 && mth4He>-2*mth3a[]+120","col"); //12C hit

  //  tree->Draw("mK3a[]+mK3a[(Iteration$+1)%3]+mK3a[(Iteration$+2)%3]:mK3a[]+mK3a[(Iteration$+1)%3]+mK3a[(Iteration$+2)%3]+mK4He>>h(100,14,18,100,0,7)","m3achF[(Iteration$+0)%3]==m3achF[(Iteration$+1)%3] && m3achR[(Iteration$+0)%3]==m3achR[(Iteration$+1)%3] && m3achF[(Iteration$+1)%3]==m3achF[(Iteration$+2)%3] && m3achR[(Iteration$+1)%3]==m3achR[(Iteration$+2)%3] && idp4He && idp3a[(Iteration$+0)%3] && idp3a[(Iteration$+1)%3] && idp3a[(Iteration$+2)%3]","colz"); //alpha 3hit

  //  tree->Draw("mK3a[]+mK3a[(Iteration$+1)%3]:mK3a[]+mK3a[(Iteration$+1)%3]+mK4He>>h(100,14,18,100,0,5)","m3achF[(Iteration$+0)%3]==m3achF[(Iteration$+1)%3] && m3achR[(Iteration$+0)%3]==m3achR[(Iteration$+1)%3] && idp4He && idp3a[(Iteration$+0)%3] && idp3a[(Iteration$+1)%3] && m3aseg[(Iteration$+2)%3]!=((m4Heseg+2)%4)","col"); //alpha 2hit

  //  tree->Draw("mK3a[]:mK3a[]+mK4He>>hh(100,14,18,100,0,5)"," idp4He && idp3a[] && m3aseg[(Iteration$+1)%3]!=((m4Heseg+2)%4) && m3aseg[(Iteration$+2)%3]!=((m4Heseg+2)%4)","col"); //alpha 1hit

  //  tree->Draw("mK12C:mK12C+mK4He>>hh(100,14,18,100,0,5)"," idp4He && idp12C","col"); //12C hit
}
