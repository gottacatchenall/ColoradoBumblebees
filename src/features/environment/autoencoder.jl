@Base.kwargs struct EnvironmentAutoencoder <: Environment
    species_encoder_outdims = 16
    shared_decoder_outdims = 128
end
