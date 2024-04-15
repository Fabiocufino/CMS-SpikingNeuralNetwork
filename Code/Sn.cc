// Loop on events for adjust thr ----------------------------------------------
    bool insert = true;

    int t_event = 0;
    int t_iev_thisepoch = 0;

    int t_EV = N_events/ 10;

    do
    {
        t_iev_thisepoch++;

        if (doprogress)
        {
            if (t_iev_thisepoch % block == 0)
            {
                cout << progress[currchar] << flush;
                currchar++;
            }
        }

        if (t_event % NROOT == 0)
        {
            last_row_event_IT = 0;
            last_row_event_OT = 0;
        }

        ReadFromProcessed(IT, OT, t_event % NROOT);

        // See if we find with track with positive latency by at least one neuron
        for (int in = 0; in < snn_in.N_neurons; in++)
        {
            Seen[pclass][in] = false; // Becomes true if the neuron in has fired for the class pclass
        }
        doneL0[pclass] = false;
        doneL1[pclass] = false;
        not_fired_bgr = true;

        // Encode hits in spike streams
        // Here we encode the position of hits through the timing of a spike,
        // In the future we might think at how the analog charge readout in the hits could also be added
        double previous_firetime = 0;
        PreSpike_Time.clear();
        PreSpike_Stream.clear();
        PreSpike_Signal.clear();

        double t_in = t_event * (max_angle + Empty_buffer) / omega; // Initial time -> every event adds 25 ns
        Encode(t_in);

        // Keep track of latency calc for each neuron in this event
        bool not_filled[snn_in.N_neurons];
        for (int in = 0; in < snn_in.N_neurons; in++)
        {
            not_filled[in] = true;
        }

        // Loop on spikes and modify neuron and synapse potentials
        // -------------------------------------------------------
        for (int ispike = 0; ispike < PreSpike_Time.size(); ispike++)
        {
            // By looping to size(), we can insert along the way and still make it to the end
            double t = PreSpike_Time[ispike];

            // Modify neuron potentials based on synapse weights
            // -------------------------------------------------
            double min_fire_time = snn_in.largenumber; // if no fire, neuron_firetime returns largenumber
            int in_first = -1;

            // Loop on neurons, but not in order to not favor any neuron
            // ---------------------------------------------------------

            // Shuffle order
            auto rng = default_random_engine{};
            shuffle(neurons_index.begin(), neurons_index.end(), rng);

            for (auto in : neurons_index)
            {
                // Compute future fire times of neurons and their order
                double fire_time = snn_in.Neuron_firetime(in, t);

                if (fire_time < min_fire_time)
                {
                    in_first = in;
                    min_fire_time = fire_time;
                }
            }
            if (in_first == -1)
                insert = true;

            // Ok, neuron in_first is going to fire next.
            // Peek at next event in list, to see if it comes before in_first fires
            // --------------------------------------------------------------------
            else
            {
                double latency = 0.;
                N_fires[in_first]++;
                snn_in.Fire_time[in_first].push_back(min_fire_time);

                // Learn weights with spike-time-dependent plasticity: long-term synaptic potentiation
                snn_in.LTP(in_first, min_fire_time, nearest_spike_approx, snn_old);

                // Reset history of this neuron
                snn_in.History_time[in_first].clear();
                snn_in.History_type[in_first].clear();
                snn_in.History_ID[in_first].clear();
                snn_in.History_time[in_first].push_back(min_fire_time);
                snn_in.History_type[in_first].push_back(0);
                snn_in.History_ID[in_first].push_back(0); // ID is not used for type 0 history events

                // IPSP for all others at relevant layer
                for (int in2 = 0; in2 < snn_in.N_neurons; in2++)
                {
                    if (in2 != in_first)
                    {
                        if (snn_in.Neuron_layer[in2] == snn_in.Neuron_layer[in_first])
                        { // inhibitions within layer or across
                            snn_in.History_time[in2].push_back(min_fire_time);
                            snn_in.History_type[in2].push_back(2);
                            snn_in.History_ID[in2].push_back(snn_in.N_InputStreams + in_first);
                        }
                    }
                }

                // Create EPS signal in L0 neuron-originated streams
                if (snn_in.Neuron_layer[in_first] == 0)
                { // this is a Layer-0 neuron
                    for (int in = snn_in.N_neuronsL[0]; in < snn_in.N_neurons; in++)
                    {
                        int is = snn_in.N_InputStreams + in_first;
                        if(!snn_in.Void_weight[in][is]){
                            snn_in.History_time[in].push_back(min_fire_time + snn_in.Delay[in][is]);
                            snn_in.History_type[in].push_back(1);
                            snn_in.History_ID[in].push_back(is);
                            snn_in.LTD(in, is, min_fire_time+ snn_in.Delay[in][in_first], nearest_spike_approx, snn_old);
                        }
                    }
                }

                // Fill spikes train histogram
                if (t_event >= N_events - 500.)
                {
                    int is = (t_event - N_events + 500) / 50;
                    double time = min_fire_time - (max_angle + Empty_buffer) / omega * (t_event / 50) * 50;
                    StreamsN[is]->Fill(time, in_first + 1);
                    if (N_part > 0)
                    {
                        fout << t_event << ", " << 2 << ", " << in_first << "," << time << "," << pclass << endl;
                    }
                }

                // Fill latency histogram
                if (N_part > 0)
                {
                    // How long did it take the first neuron to fire with respect to the arrival time of the first hit?
                    latency = min_fire_time - t_in - First_angle / omega;
                    if (latency >= 0. && not_filled[in_first])
                    {
                        if (iepoch == N_epochs - 1)
                            Latency[in_first * N_classes + pclass]->Fill(0.5 + t_iev_thisepoch, latency);
                        Seen[pclass][in_first] = true;
                        not_filled[in_first] = false;
                    }
                }
                else
                {
                    if (not_filled[in_first] && t_iev_thisepoch > NevPerEpoch * 0.9)
                    {
                        random_fire[in_first]++;
                        not_filled[in_first] = false;
                    }
                    if (in_first >= snn_in.N_neuronsL[0] && t_iev_thisepoch > NevPerEpoch * 0.9)
                    { // for Q-value calculations
                        if (not_fired_bgr)
                        {
                            atleastonefired++;
                            not_fired_bgr = false;
                        }
                    }
                }

                ispike -= 1;
                insert = false;

                // take a step back and search for another activation
                previous_firetime = min_fire_time;

            } // end if in_first fires
            // insert the new spike for the next iteration
            if (insert)
            {
                // Save information on hit-based streams for last 500 events to histograms
                if (t_event >= N_events - 500.)
                {
                    // dividing N_events in 10 groups
                    int is = (t_event - N_events + 500) / 50;
                    // time = tin + thit - tin(First event of the group)
                    double time = PreSpike_Time[ispike] - (max_angle + Empty_buffer) / omega * (t_event / 50) * 50;

                    // Histograms
                    if (PreSpike_Signal[ispike] == 1)
                    {
                        StreamsS[is]->Fill(time, PreSpike_Stream[ispike] + 1);
                        fout << t_event << ", " << PreSpike_Signal[ispike] << ", " << PreSpike_Stream[ispike] + 1 << "," << time << "," << pclass << endl;
                    }
                    else if (PreSpike_Signal[ispike] == 0)
                    {
                        StreamsB[is]->Fill(time, PreSpike_Stream[ispike] + 1);
                    }
                }
                int is = PreSpike_Stream[ispike];
                for (auto in : neurons_index)
                {
                    //  We implement a scheme where input streams produce an IE signal into L0, an EPS into L1, and L0 neurons EPS into L1
                    //  Add to neuron history, masking out L1 spikes for L0 neurons
                    if (!snn_in.Void_weight[in][is])
                    { // otherwise stream "is" does not lead to neuron "in"
                        snn_in.History_time[in].push_back(t+ snn_in.Delay[in][is]);
                        // All input spikes lead to EPSP
                        snn_in.History_type[in].push_back(1);
                        snn_in.History_ID[in].push_back(is);

                        snn_in.LTD(in, is, t+ snn_in.Delay[in][is], nearest_spike_approx, snn_old);
                    }
                }
            }
        } // end ispike loop, ready to start over

        }




        // prespike_time.push_back(time) -> time associated to an hit or to a spike coming from L0
        // prespike_Stream -> stream id associated to the hit bin or to the L0 neuron
        // prespike_Signal -> spike type: 0 if BKG, 1 if Track, 2 if NeuronFire

        // Fill efficiency histograms every NevPerEpoch events, compute Q value and Selectivity, modify parameters
        // ---------------------------------------------------------------------------------------------------

        if (t_iev_thisepoch == NevPerEpoch)
        { // we did NevPerEpoch events

            // Reset counter that inhibits efficiency and Q calculations until we reach steady state with weights
            t_iev_thisepoch = 0;
            iepoch++;
            // End of progress bar
            if (doprogress)
                cout << progress[51] << endl;

            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                for (int ic = 0; ic < N_classes; ic++)
                {
                    int combind = ic + N_classes * in;
                    Eff[combind] = fired_sum[ic][in];
                    if (gen_sum[ic] > 0)
                        Eff[combind] /= gen_sum[ic];
                    Efficiency[combind]->SetBinContent(iepoch, Eff[combind]);
                }
                float fakerate = random_fire[in] * 2. / NevPerEpoch; // there are NevPerEpoch/2 events with no tracks, where we compute random_fire per neuron
                FakeRate[in]->SetBinContent(iepoch, fakerate);
            }
            float Efftot[N_classes];
            for (int ic = 0; ic < N_classes; ic++)
            {
                float etl0 = fired_anyL0[ic];
                if (gen_sum[ic] > 0)
                    etl0 /= gen_sum[ic];
                Eff_totL0[ic]->SetBinContent(iepoch, etl0);
                float etl1 = fired_anyL1[ic]; // L1 efficiency is what counts.
                if (gen_sum[ic] > 0)
                    etl1 /= gen_sum[ic];
                Eff_totL1[ic]->SetBinContent(iepoch, etl1);
                Efftot[ic] = etl1;
            }

            selectivityL0 = Compute_Selectivity(0, 2, snn_in);
            SelectivityL0->Fill(iepoch, selectivityL0);
            selectivityL1 = Compute_Selectivity(1, 2, snn_in);
            SelectivityL1->Fill(iepoch, selectivityL1);

            // Q value is average efficiency divided by sqrt (aver eff plus aver acceptance)
            // -----------------------------------------------------------------------------
            averacctotL1 = atleastonefired * (2. / NevPerEpoch * 10.); // total acceptance, computed with 0.1*NevPerEpoch/2 events with no tracks
            averefftotL1 = 0.;
            for (int ic = 0; ic < N_classes; ic++)
            {
                averefftotL1 += Efftot[ic];
            }
            averefftotL1 /= N_classes;

            Q = Compute_Q(averefftotL1, averacctotL1, selectivityL1);

            // Fix maximum excursion of parameters with a schedule
            LR = LR_Scheduler(MaxFactor, iepoch, N_epochs);
            for (int i = 0; i < 9; i++)
            {
                max_dx[i] = LR;
            }
            for (int id = 0; id < snn_in.N_neurons * snn_in.N_streams; id++)
            {
                max_dxD[id] = 0.1 * LR;
            }

            // Re-initialize neurons
            snn_in.Init_neurons();
            // Reset hits
            Reset_hits();
            // Reset weights to initial conditions before new investigation
            snn_in.Reset_weights();
            // Init delays
            if (!updateDelays && !ReadPars && !learnDelays)
                snn_in.Init_delays(); // This unlike void connections, because we can opt to learn these at each cycle too

            cout << "         Ev. # " << t_event + 1 << " - LR = " << LR << "; Selectivity L0 = " << selectivityL0 << " L1 = " << selectivityL1
                 << "; Eff = " << averefftotL1 << " Acc = " << averacctotL1 << "; Firings: ";

            for (int in = 0; in < snn_in.N_neurons; in++)
            {
                cout << N_fires[in] << " ";
            }
            cout << endl;
            
            

            if (t_event < N_events - 1)
            { // Otherwise we graciously exit loop


                //ADJUST THRESHOLDS
                //TODO: Ema check
                // if the numer of events is 1/10 of the total number, print the firing rate of the neurons

                cout<< "\n" << endl;
                double mean_firing_rate = 0;
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    mean_firing_rate += N_fires[in];
                }
                mean_firing_rate /= snn_in.N_neurons;
                cout << "Mean firing rate is " << mean_firing_rate << "   " << endl;

                cout<<endl;

                if (update9)
                {
                    cout << "Thr Fabio -----------------------" << endl;
                    cout << " - Try TL0 = " << snn_in.Threshold[0]
                         << " TL1 = " << snn_in.Threshold[1] << " a = " << snn_in.alpha << " L1inh = " << snn_in.L1inhibitfactor
                         << " K = " << snn_in.K << " K1 = " << snn_in.K1 << " K2 = " << snn_in.K2 << " IEPC = " << snn_in.IE_Pot_const << " IPSPdf = " << snn_in.IPSP_dt_dilation << endl
                         << endl;
                }
                else
                {
                    cout << endl
                         << endl;
                }

                // Reset a few counters
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    N_fires[in] = 0.;
                }
                for (int in = 0; in < snn_in.N_neurons; in++)
                {
                    random_fire[in] = 0;
                    for (int ic = 0; ic < N_classes; ic++)
                    {
                        fired_sum[ic][in] = 0;
                    }
                }
                for (int ic = 0; ic < N_classes; ic++)
                {
                    gen_sum[ic] = 0;
                    fired_anyL0[ic] = 0;
                    fired_anyL1[ic] = 0;
                }
                atleastonefired = 0;

                // Reset progress bar
                if (doprogress)
                {
                    cout << "         " << progress[0];
                    currchar = 1;
                }
            }

        }
        t_event++; // only go to next event if we did a backward pass too
    } while (t_event < t_EV);
